from Bio import SeqIO
import sys
import os
import gzip
import numpy as np
import pandas as pd

sys.path.insert(0, '../TeloBP/')
from TeloBP import *
import constants as c

import multiprocessing as mp

def getSampleKeyFromFilename(file):
    NBtitle = [x for x in file.split('.') if "NB" in x][0]
    NBtitle = [x for x in NBtitle.split('_') if "NB" in x][0]
    NBtitle = NBtitle.replace("uq", "")
    Ftitle = [x for x in file.split('.') if "F" in x][0]
    Ftitle = [x for x in Ftitle.split('_') if "F" in x][0]
    sampleKey = Ftitle + "_" + NBtitle
    return sampleKey

# Here we read in the tables which contain read qnames and alignment information

# sampleQnamesDir = "../../Data/sampleQnamesData/assignment/"
# sampleQnamesDir = "C:/Users/Ramin Kahidi/Bioinformatics/Telomere Analysis/Nanopore project/Data/KarTongData/assignment/"
sampleQnamesDir = "./Data/KarTongData/assignment/"

sampleQnames = {}

for root, dirs, files in os.walk(sampleQnamesDir):
    for file in files:
        if not file.endswith('.txt'):
            continue
        # skip empty files
        if os.stat(os.path.join(root, file)).st_size == 0:
            continue
        print(file)
        
        sampleKey = getSampleKeyFromFilename(file)
        print(sampleKey)
        # read in table and set columns
        sampleQnames[sampleKey] = pd.read_csv(os.path.join(root, file), delimiter='\t', header=None)
        sampleQnames[sampleKey].columns = ["qname", "seqLength", "sampleQnamesChr", "subTeloAlignLength"]

# This section will read in fastq files and extract the sequences for each read

# fastqReadDir = "../../Data/Final.demultip.tagged.fastq"
fastqReadDir = "./Data/Final.demultip.tagged.fastq"

for root, dirs, files in os.walk(fastqReadDir):
    for filename in files:
        print(filename)
        if not filename.endswith(".gz") and not filename.endswith(".fastq") or "AG" in filename:
            continue

        sampleKey = getSampleKeyFromFilename(filename)
        print(sampleKey)

        if sampleKey not in sampleQnames.keys():
            continue
        print(f"Processing {filename}")
        sampleDf = sampleQnames[sampleKey]

        qnameTeloValues = []
        file = os.path.join(root, filename)
        print(file)
        if filename.endswith("fastq.gz"):
            with gzip.open(file,"rt") as handle:
                records = SeqIO.parse(handle,"fastq")
                for record in records:
                    if record.id not in sampleDf["qname"].tolist():
                        continue
                    qnameTeloValues.append([record.id, record.seq])
        elif filename.endswith("fastq"):
            for record in SeqIO.parse(file,"fastq"):
                if record.id not in sampleDf["qname"].tolist():
                    continue
                qnameTeloValues.append([record.id, record.seq])
        qnameTeloValuesDf = pd.DataFrame(qnameTeloValues, columns = ["qname", "seq"])
        sampleQnames[sampleKey] = pd.merge(sampleDf, qnameTeloValuesDf, on='qname', how='left')
        # ********************************************************************************************************************
        # Uncomment the following line when testing to break after the first file
        # break
        # ********************************************************************************************************************

# For each table in sampleQnames, remove rows with NaN in seq column
# This is because some reads may not have been present in the fastq files
popKeys = []
sampleQnamesNan = {}
for sampleKey in sampleQnames.keys():
    # print(sampleKey)
    sampleDf = sampleQnames[sampleKey]
    if "seq" not in sampleDf.keys():
        popKeys.append(sampleKey)
        continue
    # print(len(sampleDf))
    for index, row in sampleDf.iterrows():
        if row["seq"] is np.nan:
            if sampleKey not in sampleQnamesNan.keys():
                sampleQnamesNan[sampleKey] = [row]
            else:
                sampleQnamesNan[sampleKey].append(row)
            sampleDf.drop(index, inplace=True)

print(popKeys)
for key in popKeys:
    sampleQnames.pop(key)


# Find any duplicated qnames in the tables:
# No output from this cell is good. It means there are no duplicated qnames.
for sampleKey in sampleQnames.keys():
    sampleDf = sampleQnames[sampleKey]
    # print out any duplicates in the qname column
    dupQnames = sampleDf[sampleDf.duplicated(['qname'])]["qname"].tolist()
    if len(dupQnames) > 0:
        # sort by qname
        print(sampleKey)
        sortedDf = sampleDf.sort_values(by=['qname'])
        print(sortedDf)

# from pandarallel import pandarallel

outputDir = "output/TeloNPInprogress"

def rowToTeloBP(row):
    import numpy as np
    from TeloBP import getTeloNPBoundary
    # This if statement is to catch any rows which have NaN in the seq column. 
    # Ideally this should not be necessary, but it is here just in case.
    if row["seq"] is np.nan:
        return -1000
    
    teloLength = getTeloNPBoundary(row["seq"])
    return teloLength

# # The following will multiprocess the rowToTeloBP function
# for sampleKey in sampleQnames.keys():
#     sampleDf = sampleQnames[sampleKey]
    
#     pandarallel.initialize(progress_bar=True )
#     sampleDf["teloBPLengths"] = sampleDf.parallel_apply(rowToTeloBP,axis=1)

#     # I highly recommend multi-processing, but if you want to single process,
#     # comment out the above two lines and uncomment the following line:
#     # sampleDf["teloBPLengths"] = sampleDf.apply(lambda row: rowToTeloBP(row), axis=1)

#     sampleQnames[sampleKey] = sampleDf

#     # save the output to a csv file. 
#     # Note that the seq column is removed from the table before saving 
#     sampleQnames[sampleKey] = sampleQnames[sampleKey].drop(columns=["seq"])
#     sampleQnames[sampleKey].to_csv(f"{outputDir}/{sampleKey}.csv")

# The following will multiprocess each table into its own process and save the output to a csv file


def process_table(sampleKey):
    sampleDf = sampleQnames[sampleKey]
    
    # Apply the function to each row
    sampleDf["teloBPLengths"] = sampleDf.apply(rowToTeloBP, axis=1)

    # Save the output to a csv file. 
    # Note that the seq column is removed from the table before saving 
    sampleDf = sampleDf.drop(columns=["seq"])
    sampleDf.to_csv(f"{outputDir}/{sampleKey}.csv")

print(f"******************** Multiprocessing {len(sampleQnames.keys())} tables ********************")

print(f"Number of CPUs: {mp.cpu_count()}")

# Create a pool of processes
with mp.Pool(mp.cpu_count()) as pool:
    pool.map(process_table, sampleQnames.keys())

