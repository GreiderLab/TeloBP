# This script is used to trim the reference genome, removing the telomeric regions.
# Example input: python .\trimGenome.py "../Data/ncbi_dataset/data/GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna" "../outputs/trimmedGenome/GCA_009914755.4_T2T-CHM13v2.0_genomic.NoTelo.fna"

import sys

sys.path.insert(0, '..\TeloBP')

from TeloBP import trimTeloReferenceGenome

# Check if the correct number of arguments were passed
if len(sys.argv) != 3:
    print("Usage: python trimGenome.py <reference genome> <output file>")
    sys.exit()

# Call the function to trim the reference genome
# NOTE: This algorithm assumes that each chromosome is its own read, and not split into
# q and p arms.
trimTeloReferenceGenome(sys.argv[1], sys.argv[2])
