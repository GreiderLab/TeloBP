from Bio import SeqIO
import sys

from TeloBP import refRecordTeloLengths
from TeloBP.teloBoundaryHelpers import write_bed_file, recordBedData

# Check if the correct number of arguments were passed
if len(sys.argv) != 3:
    print("Usage: python teloBPBedGenome.py <reference genome> <output bed file>")
    sys.exit()
filename = sys.argv[1]
outputFile = sys.argv[2]

bed_data = []

for record in SeqIO.parse(filename, "fasta"):
    print(record.id)
    startTeloLength, endTeloLength = refRecordTeloLengths(record)
    
    recordBedData(bed_data, record, startTeloLength, endTeloLength)

write_bed_file(outputFile, bed_data)