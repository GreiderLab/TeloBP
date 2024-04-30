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
from constants import errorReturns

teloNP = False
outputColName = "teloBPLengths"

def rowToTeloBP(row, teloNP):
    import numpy as np
    from TeloBP import getTeloNPBoundary
    # This if statement is to catch any rows which have NaN in the seq column. 
    # Ideally this should not be necessary, but it is here just in case.
    if row["seq"] is np.nan:
        return -1000
    
    teloLength = -1
    if teloNP:
        teloLength = getTeloNPBoundary(row["seq"])
    else:
        teloLength = getTeloBoundary(row["seq"])
    return teloLength

def process_data(df):
    global teloNP
    df[outputColName] = df.apply(lambda row: rowToTeloBP(row, teloNP), axis=1)
    return df

# def run_analysis(dataDir, fileMode, teloNP, outputDir, progressLabel, output_frame):
def run_analysis(dataDir, fileMode, teloNPIn, outputDir):
    global teloNP 
    teloNP = teloNPIn
    print(f'Running analysis on {dataDir}')
    # Process the selected file or folder
    print(f"Data type: {'File mode' if fileMode else 'Folder mode'}")

    print(f"Use teloNP: {teloNP}")

    print("beginnging data loading")

    filenames = []
    sampleQnames = {}

    if not fileMode:
        for root, dirs, files in os.walk(dataDir):
            for filename in files:
                filenames.append((root, filename))
    else:
        filenames.append(dataDir)
        print(filenames)

    # progressLabel.config(text=f"Loading fastq files...")
    print(filenames)
    
    totalFiles = len(filenames)
    count = 0
    for filename in filenames:
        if not filename.endswith(".gz") and not filename.endswith(".fastq") or "AG" in filename:
            continue

        sampleKey = filename.split("/")[-1].split(".")[0].replace(" ", "_")
        print(f"sampling key: {sampleKey}")

        '''
        if sampleKey not in sampleQnames.keys():
            continue
        sampleDf = sampleQnames[sampleKey]
        '''

        qnameTeloValues = []
        #file = os.path.join(root, filename)
        file = filename
        print(file)
        if filename.endswith("fastq.gz"):
            with gzip.open(file,"rt") as handle:
                records = SeqIO.parse(handle,"fastq")
                for record in records:
                    '''
                    if record.id not in sampleDf["qname"].tolist():
                        continue
                    '''
                    qnameTeloValues.append([record.id, record.seq])
        elif filename.endswith("fastq"):
            for record in SeqIO.parse(file,"fastq"):
                '''
                if record.id not in sampleDf["qname"].tolist():
                    continue
                '''

                qnameTeloValues.append([record.id, record.seq])
        qnameTeloValuesDf = pd.DataFrame(qnameTeloValues, columns = ["qname", "seq"])
        #sampleQnames[sampleKey] = pd.merge(sampleDf, qnameTeloValuesDf, on='qname', how='left')
        sampleQnames[sampleKey] = qnameTeloValuesDf
        
        count += 1
        # progressLabel.config(text=f"Loading fastq files... ({count}/{totalFiles})")

    print("data loaded")
    # progressLabel.config(text=f"data loaded")


    #global outputDir
    # if totalFiles <= 10:
    if fileMode:
        print("RUNNING SINGLE FILE")
        
        count = 0
        # The following will multiprocess the rowToTeloBP function
        for sampleKey in sampleQnames.keys():
            sampleDf = sampleQnames[sampleKey]
            inputLen = len(sampleDf)

            # def rowToTeloBP(row, teloNP):
            #     import numpy as np
            #     from TeloBP import getTeloNPBoundary
            #     # This if statement is to catch any rows which have NaN in the seq column. 
            #     # Ideally this should not be necessary, but it is here just in case.
            #     if row["seq"] is np.nan:
            #         return -1000
                
            #     teloLength = -1
            #     if teloNP:
            #         teloLength = getTeloNPBoundary(row["seq"])
            #     else:
            #         teloLength = getTeloBoundary(row["seq"])
            #     return teloLength

            # pandarallel.initialize(progress_bar=True)
            # sampleDf[outputColName] = sampleDf.parallel_apply(lambda row: rowToTeloBP(row, teloNP), axis=1)


            # Split DataFrame into chunks
            n_chunks = mp.cpu_count()  # Number of chunks is equal to the number of cores
            df_chunks = np.array_split(sampleDf, n_chunks)

            # Create a multiprocessing Pool
            pool = mp.Pool(n_chunks)

            # Process each chunk in parallel
            df_chunks = pool.map(process_data, df_chunks)

            # Concatenate the chunks back into a single DataFrame
            sampleDf = pd.concat(df_chunks)

            # I highly recommend multi-processing, but if you want to single process,
            # comment out the above two lines and uncomment the following line:
            # sampleDf[outputColName] = sampleDf.apply(lambda row: rowToTeloBP(row), axis=1)

            sampleQnames[sampleKey] = sampleDf

            print(f"inputLen: {inputLen}")
            print(f"outputLen: {len(sampleDf)}")
            if inputLen != len(sampleDf):
                print(f"Error: input length {inputLen} does not match output length {len(sampleDf)}")

            # save the output to a csv file. 
            # Note that the seq column is removed from the table before saving 
            sampleQnames[sampleKey] = sampleQnames[sampleKey].drop(columns=["seq"])
            dfNoNeg = sampleQnames[sampleKey][sampleQnames[sampleKey][outputColName] > 0]
            dfNoNeg.reset_index(drop=True, inplace=True)
            dfNoNeg.to_csv(f"{outputDir}/{sampleKey}.csv")
            count += 1
            # progressLabel.config(text=f"Processing tables... ({count}/{totalFiles})")
    else:
        print("RUNNING MULTI-PROCESSING")
        count = 0
        def process_table(sampleKey):
            sampleDf = sampleQnames[sampleKey]
            
            # Apply the function to each row
            sampleDf[outputColName] = sampleDf.apply(lambda row: rowToTeloBP(row, teloNP), axis=1)
            # sampleDf[outputColName] = sampleDf.parallel_apply(lambda row: rowToTeloBP(row, teloNP), axis=1)

            # Save the output to a csv file. 
            # Note that the seq column is removed from the table before saving 
            sampleDf = sampleDf.drop(columns=["seq"])
            dfNoNeg = sampleQnames[sampleKey][sampleQnames[sampleKey][outputColName] > 0]
            dfNoNeg.reset_index(drop=True, inplace=True)
            dfNoNeg.to_csv(f"{outputDir}/{sampleKey}.csv")
            nonlocal count
            count += 1
            # progressLabel.config(text=f"Processing tables... ({count}/{totalFiles})")

        print(f"******************** Multiprocessing {len(sampleQnames.keys())} tables ********************")
        print(f"Number of CPUs: {mp.cpu_count()}")

        # Create a pool of processes
        with mp.Pool(mp.cpu_count()) as pool:
            pool.map(process_table, sampleQnames.keys())

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
        print(sampleKey)
        print(len(sampleDf[sampleDf[outputColName] < 0]))
        print(f"Total reads: {len(sampleDf)}")
        print(f"Percentage of reads with teloBP < 0: {len(sampleDf[sampleDf['teloBPLengths'] < 0])/len(sampleDf)*100}%")
        print(f"Average teloBP length: {sampleDf['teloBPLengths'].mean()}")
        print(f"Average teloBP length (excluding -ve values): {sampleDf[sampleDf['teloBPLengths'] > 0]['teloBPLengths'].mean()}")

        totalReads += len(sampleDf)
        totalReadLengths += sampleDf[outputColName].sum()
        nonNegativeReads += len(sampleDf[sampleDf[outputColName] > 0])
        initErrors += len(sampleDf[sampleDf[outputColName] == errorReturns["init"]])
        fusedReadErrors += len(sampleDf[sampleDf[outputColName] == errorReturns["fusedRead"]])
        strandTypeErrors += len(sampleDf[sampleDf[outputColName] == errorReturns["strandType"]])
        seqNotFound += len(sampleDf[sampleDf[outputColName] == errorReturns["seqNotFound"]])

    avgTeloBP = totalReadLengths / nonNegativeReads

    if totalReads == (nonNegativeReads + initErrors + fusedReadErrors + strandTypeErrors + seqNotFound):
        print("Total reads match")
    else: 
        print(f"Total reads: {totalReads}")
        print(f"Non-negative reads: {nonNegativeReads}")
        print(f"Initialization errors: {initErrors}")
        print(f"Fused read errors: {fusedReadErrors}")
        print(f"Strand type errors: {strandTypeErrors}")
        print(f"Seq not found errors: {seqNotFound}")

        raise Warning("Total reads do not match!!!!!")

    # Data to plot
    #labels = 'Python', 'C++', 'Ruby', 'Java'
    labels = 'Passed Reads', 'initialization error', 'fused read error', 'strand type error', 'seq Not Found'
    sizes = [nonNegativeReads, initErrors, fusedReadErrors, strandTypeErrors, seqNotFound]
    colors = ['green', 'orange', 'yellow', 'blue','red']

    outputText = f"Total reads: {totalReads}\n" + \
        f"Average telo length: {avgTeloBP}\n" + \
        f"Non-negative reads: {nonNegativeReads}\n" + \
        f"Initialization errors: {initErrors}\n" + \
        f"Fused read errors: {fusedReadErrors}\n" + \
        f"Strand type errors: {strandTypeErrors}\n" + \
        f"Seq not found errors: {seqNotFound}\n" 

    # graph(labels, sizes, colors, outputText, output_frame)



if __name__ == "__main__":
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Command line input handling for run_analysis function')

    # Add the arguments
    parser.add_argument('dataDir', type=str, help='Path to the data directory')
    parser.add_argument('outputDir', type=str, help='Path to the output directory')
    parser.add_argument('--fileMode', action='store_true', help='Flag to indicate if we are looking at a single file or a direcotry of files')
    parser.add_argument('--teloNPBool', action='store_true', help='Flag to indicate whether to use teloNP')
    # parser.add_argument('--progressLabel', type=str, help='Progress label')
    # parser.add_argument('--output_frame', type=str, help='Output frame')

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the run_analysis function with the parsed arguments
    run_analysis(args.dataDir, args.fileMode, args.teloNPBool, args.outputDir)
    # run_analysis("../data", False, True, "../output", None, None)