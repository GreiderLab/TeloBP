from TeloBP import trimTeloReferenceGenome
import sys

sys.path.insert(0, '../TeloBP')

print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))

trimTeloReferenceGenome(sys.argv[1], sys.argv[2])