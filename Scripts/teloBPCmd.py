# This script is used to run the teloBP analysis on a set of fastq files.
# The script can be run from the command line using the following command:
# python3 teloBPCmd.py <dataDir> <outputDir> [--fileMode] [--teloNP] [--chemistry R10|R9|R9_archive] [--workers N] [-v]
# The script takes the following arguments:
# dataDir: The path to the data directory containing the fastq files to be analyzed. 
#          OR the path to a single fastq file to be analyzed. set the --fileMode flag to indicate this.
# outputDir: The path to the output directory where the results will be saved.
# --fileMode: Flag indicating that a single file is being analyzed.
# --teloNP: Flag indicating that the teloNP analysis should be run instead of the teloBP analysis.
# --chemistry: The nanopore chemistry TeloNP settings are chosen for. R10 (default), R9 for Guppy-basecalled R9.4.1 reads,
#              or R9_archive to reproduce results from TeloNP v1.
# --targetQnamesCSV: Path to a csv file containing the qnames of the reads to be analyzed. The analysis will only be run on these reads.
# --save_graphs: Flag indicating to generate graphs during analysis and save them to a pdf file. They will be saved in the main directory.
# The program will run in single threaded mode when --save_graphs is set to True.
# --workers: The number of worker processes the reads are analyzed with. Defaults to all CPUs.
# -v: Flag to enable verbose output.


import os
from Bio import SeqIO
import gzip
import numpy as np
import pandas as pd
import sys
import argparse
import json
import platform
import subprocess
import time
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import Bio
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import TeloBP as teloBPPackage

from TeloBP import *
# from constants import errorReturns
errorReturns = {"init": -1, "fusedRead": -10,"strandType": -20, "seqNotFound": -1000}

# Note: with the spawn start method (Windows, macOS) every worker process imports this script again,
# so module level code must not have side effects like creating files.
teloNP = False
chemistryGlobal = "R10"
# The graphs pdf is only opened by run_analysis when graphs are saved, which always runs in a single process
pdf = None
showGraphsGlobal = False
outputColName = "teloBPLengths"
workersGlobal = 1
# Seconds between "processed X/Y reads" progress lines
progressInterval = 5

verbose = False

def vprint(*args, **kwargs):
    if verbose:
        print(*args, **kwargs)

# def rowToTeloBP(row, teloNP):
#     import numpy as np
#     from TeloBP import getTeloNPBoundary, getTeloBoundary
#     # This if statement is to catch any rows which have NaN in the seq column. 
#     # Ideally this should not be necessary, but it is here just in case.
#     if row["seq"] is np.nan:
#         return errorReturns["seqNotFound"]
    
#     teloLength = -1
#     if teloNP:
#         teloLength = getTeloNPBoundary(row["seq"])
#     else:
#         teloLength = getTeloBoundary(row["seq"])
#     return teloLength

# rowToTeloBP and rowToTeloNP run in the worker processes, so they have to stay top level functions that can be pickled.
# row is a dict with "qname" and "seq".
def rowToTeloBP(row, sampleKey, showGraphsGlobal):
    if row["seq"] is np.nan:
        return errorReturns["seqNotFound"]
    teloLength = -1
    
    pdfLocal = None
    if showGraphsGlobal:
        global pdf
        plt.figure()
        plt.text(0.5, 0.5, f'Sample Key: {sampleKey}, Qname: {row["qname"]}', ha='center', va='center', fontsize=12)
        plt.axis('off')
        pdf.savefig()  # Save the text figure
        pdfLocal = pdf
    teloLength = getTeloBoundary(row["seq"], showGraphs = showGraphsGlobal, pdf = pdfLocal)
    return teloLength

def rowToTeloNP(row, sampleKey, showGraphsGlobal, chemistry):
    if row["seq"] is np.nan:
        return errorReturns["seqNotFound"]
    teloLength = -1

    pdfLocal = None
    if showGraphsGlobal:
        plt.figure()
        plt.text(0.5, 0.5, f'Sample Key: {sampleKey}, Qname: {row["qname"]}', ha='center', va='center', fontsize=9)
        plt.axis('off')
        pdf.savefig()  # Save the text figure
        pdfLocal = pdf
    teloLength = getTeloNPBoundary(row["seq"],showGraphs = showGraphsGlobal, pdf = pdfLocal, chemistry = chemistry)
    return teloLength

def getDefaultWorkers():
    # All the CPUs this process may use. On Linux this respects CPU limits set by e.g. SLURM or taskset.
    if hasattr(os, "sched_getaffinity"):
        cpus = len(os.sched_getaffinity(0))
    else:
        cpus = os.cpu_count() or 1
    # Windows supports at most 61 worker processes
    if sys.platform == "win32":
        cpus = min(cpus, 61)
    return cpus

def getChunksize(numReads, workers):
    # Reads are sent to the workers in chunks, so thousands of short reads aren't slowed down by per-read overhead.
    # Chunks are kept small (about 8 per worker, at most 50 reads) so that the workers finish at about the same time.
    return max(1, min(50, numReads // (workers * 8)))

def printProgress(teloLengths, total, sampleKey):
    # Passes the results through, printing "processed X/Y reads" every progressInterval seconds and once at the end
    lastPrint = time.monotonic()
    for count, teloLength in enumerate(teloLengths, 1):
        yield teloLength
        if count == total or time.monotonic() - lastPrint >= progressInterval:
            print(f"{sampleKey}: processed {count}/{total} reads", flush=True)
            lastPrint = time.monotonic()

def process_sampleDf(sampleQnames, sampleKey, outputDir, executor=None):
    # executor is the pool of worker processes. When it is None the reads are analyzed in this process.
    sampleDf = sampleQnames[sampleKey]
    inputLen = len(sampleDf)
    
    # get current path
    currentPath = os.path.dirname(os.path.realpath(__file__))
    vprint(f"Current path: {currentPath}")

    if teloNP:
        rowFunction = partial(rowToTeloNP, sampleKey=sampleKey, showGraphsGlobal=showGraphsGlobal, chemistry=chemistryGlobal)
    else:
        rowFunction = partial(rowToTeloBP, sampleKey=sampleKey, showGraphsGlobal=showGraphsGlobal)
    # Plain dicts are much cheaper to send to the worker processes than pandas rows
    rows = sampleDf.to_dict("records")

    if executor is None:
        print(f"Analyzing {inputLen} reads from {sampleKey} in a single process", flush=True)
        # Graphs are always drawn here, so they all go into the one pdf
        teloLengths = map(rowFunction, rows)
    else:
        print(f"Analyzing {inputLen} reads from {sampleKey} with {workersGlobal} worker processes", flush=True)
        # map returns the results in the same order as the reads
        teloLengths = executor.map(rowFunction, rows, chunksize=getChunksize(inputLen, workersGlobal))
    sampleDf[outputColName] = list(printProgress(teloLengths, inputLen, sampleKey))

    vprint(f"inputLen: {inputLen}")
    vprint(f"outputLen: {len(sampleDf)}")
    if inputLen != len(sampleDf):
        print(f"Error: starting length of {inputLen} does not match output table length {len(sampleDf)}")

    # save the output to a csv file. 
    # Note that the seq column is removed from the table before saving 
    sampleDf = sampleDf.drop(columns=["seq"])
    
    saveDfNoNeg(sampleDf, outputColName, outputDir, sampleKey)
    return sampleDf

def process_data(df):
    global teloNP
    df[outputColName] = df.apply(lambda row: rowToTeloBP(row, teloNP), axis=1)
    return df

def getGitCommit():
    # The commit of the TeloBP checkout being run, with "-dirty" if it has uncommitted changes.
    # None when TeloBP wasn't installed from a git checkout.
    packageDir = os.path.dirname(os.path.realpath(teloBPPackage.__file__))
    try:
        commit = subprocess.run(["git", "-C", packageDir, "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", packageDir, "status", "--porcelain", "--untracked-files=no"], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None
    return commit + ("-dirty" if dirty else "")

def writeRunMetadata(outputDir, filenames, sampleQnames):
    # Saves the version and settings of this run next to its results, so they can be reproduced
    metadata = {
        "teloBPVersion": teloBPPackage.__version__,
        "gitCommit": getGitCommit(),
        "runAt": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "method": "TeloNP" if teloNP else "TeloBP",
        "settings": getTeloNPSettings(chemistryGlobal) if teloNP else getTeloBPSettings(),
        # Number of processes the reads were analyzed with (1 when graphs were saved). It doesn't change the results.
        "workers": workersGlobal,
        "inputFiles": filenames,
        "readsPerSample": {sampleKey: len(sampleDf) for sampleKey, sampleDf in sampleQnames.items()},
        "python": platform.python_version(),
        "packages": {"biopython": Bio.__version__, "numpy": np.__version__, "pandas": pd.__version__},
    }
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    with open(f"{outputDir}/teloBP_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

def saveDfNoNeg(df, outputColName, outputDir, sampleKey):
    dfNoNeg = df[df[outputColName] > 0]
    dfNoNeg = dfNoNeg.reset_index(drop=True)
    # create outdir if it doesn't exist
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    dfNoNeg.to_csv(f"{outputDir}/{sampleKey}.csv")

# def run_analysis(dataDir, fileMode, teloNP, outputDir, progressLabel, output_frame):
def run_analysis(dataDir, fileMode, teloNPIn, outputDir, save_graphs=False, targetQnamesCSV=None, chemistry="R10", workers=None):
    # workers: number of worker processes, None for all CPUs
    global teloNP, chemistryGlobal, workersGlobal, pdf
    teloNP = teloNPIn
    chemistryGlobal = chemistry
    if teloNP:
        global outputColName
        outputColName = "teloNPLengths"
    vprint(f'Running analysis on {dataDir}')
    # Process the selected file or folder
    vprint(f"Data type: {'File mode' if fileMode else 'Folder mode'}")
    print(f"Use teloNP: {teloNP}" + (f" ({chemistryGlobal} settings)" if teloNP else ""))
    vprint("Beginnging data loading")

    if targetQnamesCSV is None and save_graphs:
        print("WARNING: save_graphs flag is set but no targetQnamesCSV file is provided. This configuration would likely crash the program, so the save_graphs flag will be set to False. If this was your intention, please modify the code.")
        save_graphs = False
    elif targetQnamesCSV is not None and save_graphs:
        global showGraphsGlobal
        showGraphsGlobal = True
    
    if workers is None:
        workers = getDefaultWorkers()
    if workers < 1:
        raise ValueError(f"workers must be at least 1, got {workers}")
    if showGraphsGlobal:
        # Every graph goes into one pdf, so the reads are analyzed in this process only
        pdf = PdfPages('teloGraphs.pdf')
        workers = 1
    workersGlobal = workers
    vprint(f"Worker processes: {workersGlobal}")

    targetQnames = None
    if targetQnamesCSV is not None:
        targetQnames = pd.read_csv(targetQnamesCSV, header=None)
        targetQnames = set(targetQnames[0].to_list())
    
    filenames = []
    sampleQnames = {}

    if not fileMode:
        for root, dirs, files in os.walk(dataDir):
            for filename in files:
                if filename.endswith(".gz") or filename.endswith(".fastq"):
                    filenames.append(os.path.join(root, filename))
    else:
        filenames.append(dataDir)

    # progressLabel.config(text=f"Loading fastq files...")
    vprint(f"filenames: {filenames}")
    
    totalFiles = len(filenames)
    count = 0
    # The worker processes are started once and reused for every file
    executor = None
    if workersGlobal > 1:
        executor = ProcessPoolExecutor(max_workers=workersGlobal)
    try:
        for filename in filenames:
            if not filename.endswith(".gz") and not filename.endswith(".fastq"):
                continue

            sampleKey = "_".join(filename.split("/")[-1].split("\\")[-1].split(".")[:-1]).replace(" ", "_")
            vprint(f"sampling key: {sampleKey}")

            qnameTeloValues = []
            #file = os.path.join(root, filename)
            file = filename
            vprint(file)
            if filename.endswith("fastq.gz"):
                with gzip.open(file,"rt") as handle:
                    records = SeqIO.parse(handle,"fastq")
                    for record in records:
                        if targetQnames is not None and record.id not in targetQnames:
                            continue
                        print(record.id)
                        qnameTeloValues.append([record.id, record.seq])
            elif filename.endswith("fastq"):
                try:
                    for record in SeqIO.parse(file,"fastq"):
                        if targetQnames is not None and record.id not in targetQnames:
                            continue
                        print(record.id)
                        qnameTeloValues.append([record.id, record.seq])
                except Exception as e:
                    print(f"Error while reading file: {file}")
                    print(e)
            else: 
                continue

            qnameTeloValuesDf = pd.DataFrame(qnameTeloValues, columns = ["qname", "seq"])

            if len(qnameTeloValuesDf) == 0:
                vprint(f"No records in file: {file}")
                continue

            #sampleQnames[sampleKey] = pd.merge(sampleDf, qnameTeloValuesDf, on='qname', how='left')
            sampleQnames[sampleKey] = qnameTeloValuesDf

            sampleQnames[sampleKey] = process_sampleDf(sampleQnames, sampleKey, outputDir, executor)

            count += 1
        
            # progressLabel.config(text=f"Loading fastq files... ({count}/{totalFiles})")
    finally:
        if executor is not None:
            executor.shutdown(cancel_futures=True)

    # Save pdf
    if pdf is not None:
        pdf.close()

    writeRunMetadata(outputDir, filenames, sampleQnames)

    # ********** Stats **********

    totalReads = 0
    totalReadLengths = 0
    nonNegativeReads = 0
    initErrors = 0
    fusedReadErrors = 0
    strandTypeErrors = 0
    seqNotFound = 0


    # go through table and see how many lengths are less than 0
    for sampleKey in sampleQnames.keys():
        sampleDf = sampleQnames[sampleKey]
        sampleTotalReads = len(sampleDf)
        sampleNonNegativeReads = len(sampleDf[sampleDf[outputColName] >= 0])
        sampleInitErrors = len(sampleDf[sampleDf[outputColName] == errorReturns["init"]])
        sampleFusedReadErrors = len(sampleDf[sampleDf[outputColName] == errorReturns["fusedRead"]])
        sampleStrandTypeErrors = len(sampleDf[sampleDf[outputColName] == errorReturns["strandType"]])
        sampleSeqNotFound = len(sampleDf[sampleDf[outputColName] == errorReturns["seqNotFound"]])

        vprint(f"Sample: {sampleKey}")
        vprint(f"Average teloBP length: {sampleDf[sampleDf[outputColName] >= 0][outputColName].mean()}")

        vprint(f"\nTotal reads (including errors): {sampleTotalReads}")
        vprint(f"Reads without errors: {sampleNonNegativeReads}")
        vprint(f"Initialization errors: {sampleInitErrors}")
        vprint(f"Fused read errors: {sampleFusedReadErrors}")
        vprint(f"Strand type errors: {sampleStrandTypeErrors}") 
        vprint(f"Seq not found errors: {sampleSeqNotFound}")

        totalReads += sampleTotalReads
        totalReadLengths += sampleDf[sampleDf[outputColName] >= 0][outputColName].sum()
        nonNegativeReads += sampleNonNegativeReads
        initErrors += sampleInitErrors
        fusedReadErrors += sampleFusedReadErrors
        strandTypeErrors += sampleStrandTypeErrors
        seqNotFound += sampleSeqNotFound

    if nonNegativeReads == 0:
        avgTeloBP = 0
    else:
        avgTeloBP = totalReadLengths / nonNegativeReads

    if totalReads == (nonNegativeReads + initErrors + fusedReadErrors + strandTypeErrors + seqNotFound):
        vprint("Total reads match")
    else: 
        print(f"Total reads: {totalReads}")
        print(f"Average telo length: {avgTeloBP}")
        print(f"Non-negative reads: {nonNegativeReads}")
        print(f"Initialization errors: {initErrors}")
        print(f"Fused read errors: {fusedReadErrors}")
        print(f"Strand type errors: {strandTypeErrors}")
        print(f"Seq not found errors: {seqNotFound}")

        raise Warning("Error: Total reads do not match. Likely a script error.")

    # Data to plot
    # labels = 'Passed Reads', 'initialization error', 'fused read error', 'strand type error', 'seq Not Found'
    # sizes = [nonNegativeReads, initErrors, fusedReadErrors, strandTypeErrors, seqNotFound]
    # colors = ['green', 'orange', 'yellow', 'blue','red']

    # outputText = f"Total reads: {totalReads}\n" + \
    #     f"Average telo length: {avgTeloBP}\n" + \
    #     f"Non-negative reads: {nonNegativeReads}\n" + \
    #     f"Initialization errors: {initErrors}\n" + \
    #     f"Fused read errors: {fusedReadErrors}\n" + \
    #     f"Strand type errors: {strandTypeErrors}\n" + \
    #     f"Seq not found errors: {seqNotFound}\n" 

    # graph(labels, sizes, colors, outputText, output_frame)



if __name__ == "__main__":
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Command line input handling for run_analysis function')

    # Add the arguments
    parser.add_argument('dataDir', type=str, help='Path to the data directory')
    parser.add_argument('outputDir', type=str, help='Path to the output directory')
    parser.add_argument('--fileMode', action='store_true', help='Flag to indicate if we are looking at a single file or a direcotry of files')
    parser.add_argument('--teloNP', action='store_true', help='Flag to indicate whether to use teloNP')
    parser.add_argument('--chemistry', choices=['R10', 'R9', 'R9_archive'], default='R10', help='Nanopore chemistry for the teloNP settings: R10 (default), R9 for Guppy-basecalled R9.4.1 reads, or R9_archive to reproduce TeloNP v1 results')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--save_graphs', action='store_true', help='Flag to indicate whether to save graphs')
    parser.add_argument('--targetQnamesCSV', type=str, help='Path to an input csv that has the qnammes we want to test')  # Optional argument for input file
    parser.add_argument('--workers', type=int, default=None, help='Number of worker processes to analyze the reads with (default: all CPUs). Graphs are always made in a single process.')

    # parser.add_argument('--progressLabel', type=str, help='Progress label')

    # Parse the command line arguments
    args = parser.parse_args()
    verbose = args.verbose
    if args.workers is not None and args.workers < 1:
        parser.error("--workers must be at least 1")

    # Call the run_analysis function with the parsed arguments
    run_analysis(args.dataDir, args.fileMode, args.teloNP, args.outputDir, save_graphs=args.save_graphs, targetQnamesCSV=args.targetQnamesCSV, chemistry=args.chemistry, workers=args.workers)
    # run_analysis("../data", False, True, "../output", None, None)
