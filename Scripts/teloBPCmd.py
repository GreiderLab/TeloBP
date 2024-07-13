import os
from Bio import SeqIO
import gzip
import numpy as np
import pandas as pd
import multiprocessing as mp
from pandarallel import pandarallel
import sys
import argparse

sys.path.insert(0, '../TeloBP')
from TeloBP import *
# from constants import errorReturns
errorReturns = {"init": -1, "fusedRead": -10,"strandType": -20, "seqNotFound": -1000}

teloNP = False
outputColName = "teloBPLengths"

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

def rowToTeloBP(row):
    import sys
    sys.path.insert(0, '../TeloBP')
    import numpy as np
    from TeloBP import getTeloBoundary
    if row["seq"] is np.nan:
        return errorReturns["seqNotFound"]
    teloLength = -1
    teloLength = getTeloBoundary(row["seq"])
    return teloLength

def rowToTeloNP(row):
    import sys
    sys.path.insert(0, '../TeloBP')
    import numpy as np
    from TeloBP import getTeloNPBoundary
    if row["seq"] is np.nan:
        return errorReturns["seqNotFound"]
    teloLength = -1
    teloLength = getTeloNPBoundary(row["seq"])
    return teloLength

def process_sampleDf(sampleQnames, sampleKey, outputDir):
    sampleDf = sampleQnames[sampleKey]
    inputLen = len(sampleDf)
    
    # get current path
    currentPath = os.path.dirname(os.path.realpath(__file__))
    vprint(f"Current path: {currentPath}")

    pandarallel.initialize(progress_bar=True)
    if teloNP:
        sampleDf[outputColName] = sampleDf.parallel_apply(rowToTeloNP,axis=1)
    else:
        sampleDf[outputColName] = sampleDf.parallel_apply(rowToTeloBP,axis=1)

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

def saveDfNoNeg(df, outputColName, outputDir, sampleKey):
    dfNoNeg = df[df[outputColName] > 0]
    dfNoNeg.reset_index(drop=True, inplace=True)
    # create outdir if it doesn't exist
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    dfNoNeg.to_csv(f"{outputDir}/{sampleKey}.csv")

# def run_analysis(dataDir, fileMode, teloNP, outputDir, progressLabel, output_frame):
def run_analysis(dataDir, fileMode, teloNPIn, outputDir):
    global teloNP 
    teloNP = teloNPIn
    if teloNP:
        global outputColName
        outputColName = "teloNPLengths"
    vprint(f'Running analysis on {dataDir}')
    # Process the selected file or folder
    vprint(f"Data type: {'File mode' if fileMode else 'Folder mode'}")
    print(f"Use teloNP: {teloNP}")
    vprint("Beginnging data loading")

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
                    qnameTeloValues.append([record.id, record.seq])
        elif filename.endswith("fastq"):
            try:
                for record in SeqIO.parse(file,"fastq"):
                    qnameTeloValues.append([record.id, record.seq])
            except Exception as e:
                print(f"Error while reading file: {file}")
                print(e)
        else: 
            continue
        qnameTeloValuesDf = pd.DataFrame(qnameTeloValues, columns = ["qname", "seq"])
        #sampleQnames[sampleKey] = pd.merge(sampleDf, qnameTeloValuesDf, on='qname', how='left')
        sampleQnames[sampleKey] = qnameTeloValuesDf

        sampleQnames[sampleKey] = process_sampleDf(sampleQnames, sampleKey, outputDir)

        count += 1
        
        # progressLabel.config(text=f"Loading fastq files... ({count}/{totalFiles})")

    # The following code does not work at the moment. 
    # else:
    #     vprint("RUNNING MULTI-PROCESSING")
    #     count = 0
    #     def process_table(sampleKey):
    #         sampleDf = sampleQnames[sampleKey]
            
    #         # Apply the function to each row
    #         sampleDf[outputColName] = sampleDf.apply(lambda row: rowToTeloBP(row, teloNP), axis=1)

    #         # Note that the seq column is removed from the table before saving 
    #         sampleDf = sampleDf.drop(columns=["seq"])
    #         sampleQnames[sampleKey] = sampleDf

    #         # Save the output to a csv file. 
    #         saveDfNoNeg(sampleDf, outputColName, outputDir, sampleKey)
    #         nonlocal count
    #         count += 1
    #         # progressLabel.config(text=f"Processing tables... ({count}/{totalFiles})")

    #     vprint(f"******************** Multiprocessing {len(sampleQnames.keys())} tables ********************")
    #     vprint(f"Number of CPUs: {mp.cpu_count()}")

    #     # Create a pool of processes
    #     with mp.Pool(mp.cpu_count()) as pool:
    #         pool.map(process_table, sampleQnames.keys())


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
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    
    # parser.add_argument('--progressLabel', type=str, help='Progress label')

    # Parse the command line arguments
    args = parser.parse_args()
    verbose = args.verbose

    # Call the run_analysis function with the parsed arguments
    run_analysis(args.dataDir, args.fileMode, args.teloNP, args.outputDir)
    # run_analysis("../data", False, True, "../output", None, None)
