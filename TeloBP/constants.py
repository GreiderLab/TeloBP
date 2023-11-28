# These are constants used in the TeloBP package.

# These are the default telomere compositions for the teloBP function.
# Format: [pattern, expected % of telomere composed of pattern, length of pattern(required for Regex patterns)]
# Examples:
#     ["T", 2/6],
#     ["A", 1/6],
#     ["G", 3/6],
#     ["C", 0/6],
#     ["GGG", 3/6],
#     ["TT", 2/6],
#     ["GGG|AAA", 3/6, 3],
# expectedTeloCompositionQ = [
#     ["T", 2/6],
#     ["A", 1/6],
#     ["G", 3/6],
#     ["C", 0/6],
#     ["GGG", 3/6],
#     ["TT", 2/6],
#     ["GGG|AAA", 3/6, 3],
# ]

# expectedTeloCompositionP = [
#     ["T", expectedTeloCompositionQ[0][1]],
#     ["A", expectedTeloCompositionQ[1][1]],
#     ["G", expectedTeloCompositionQ[2][1]],
#     ["C", expectedTeloCompositionQ[3][1]],
#     ["CCC", expectedTeloCompositionQ[4][1]],
#     ["AA", expectedTeloCompositionQ[5][1]],
#     ["CCCTAA|CTTCTT|CCCTGG|CCTGG", 6/6, 6],
#     ["CCC|TTT", expectedTeloCompositionQ[6][1], 3]
# ]

expectedTeloCompositionQ = [
    ["GGG", 3/6],
]

expectedTeloCompositionP = [
    ["CCC", expectedTeloCompositionQ[0][1]],
]

# When the area under the curve of the offsets is calculated, this
# constant is used to mark the point where the slope of the curve
# starts to plateau. Basically, once the difference between two
# values in the area under the curve is less than this, the slope
# is considered to be plateauing.
areaDiffsThreshold = 0.2


# The following dictionaries are used to test the telomere length and
# position given by teloBP, and were created by manual inspection
manualLabelsCHM13Positions = {
    "chr_01p": 2704,
    "chr_01q": 248384118,
    "chr_02p": 3617,
    "chr_02q": 242694130,
    "chr_03p": 2640,
    "chr_03q": 201101333,
    "chr_04p": 3267,
    "chr_04q": 193572624,
    "chr_05p": 2295,
    "chr_05q": 182043909,
    "chr_06p": 2896,
    "chr_06q": 172123850,
    "chr_07p": 3374,
    "chr_07q": 160565206,
    "chr_08p": 2517,
    "chr_08q": 146256710,
    "chr_09p": 3586,
    "chr_09q": 150614270,
    "chr_10p": 2637,
    "chr_10q": 134754931,
    "chr_11p": 1987,
    "chr_11q": 135125180,
    "chr_12p": 3101,
    "chr_12q": 133322209,
    "chr_13p": 2543,
    "chr_13q": 113563180,
    "chr_14p": 2074,
    "chr_14q": 101159840,
    "chr_15p": 3264,
    "chr_15q": 99750266,
    "chr_16p": 2323,
    "chr_16q": 96327691,
    "chr_17p": 2209,
    "chr_17q": 84273876,
    "chr_18p": 2015,
    "chr_18q": 80539049,
    "chr_19p": 2284,
    "chr_19q": 61704424,
    "chr_20p": 2727,
    "chr_20q": 66207102,
    "chr_21p": 3012,
    "chr_21q": 45086109,
    "chr_22p": 4577,
    "chr_22q": 51321952,
    "chr_Xp": 1826,
    "chr_Xq": 154256623,
    "chr_Yp": 5658,
    "chr_Yq": 62453641,
}

manualLabelsCHM13 = {
    "chr01p": 2704,
    "chr01q": 3210,
    "chr02p": 3617,
    "chr02q": 2622,
    "chr03p": 2640,
    "chr03q": 4615,
    "chr04p": 3267,
    "chr04q": 2321,
    "chr05p": 2295,
    "chr05q": 1530,
    "chr06p": 2896,
    "chr06q": 2778,
    "chr07p": 3374,
    "chr07q": 2222,
    "chr08p": 2517,
    "chr08q": 2621,
    "chr09p": 3586,
    "chr09q": 2977,
    "chr10p": 2637,
    "chr10q": 3203,
    "chr11p": 1987,
    "chr11q": 2589,
    "chr12p": 3101,
    "chr12q": 2339,
    "chr13p": 2543,
    "chr13q": 3506,
    "chr14p": 2074,
    "chr14q": 1652,
    "chr15p": 3264,
    "chr15q": 2929,
    "chr16p": 2323,
    "chr16q": 2683,
    "chr17p": 2209,
    "chr17q": 3021,
    "chr18p": 2015,
    "chr18q": 3489,
    "chr19p": 2284,
    "chr19q": 2940,
    "chr20p": 2727,
    "chr20q": 3153,
    "chr21p": 3012,
    "chr21q": 4573,
    "chr22p": 4577,
    "chr22q": 2974,
    "chrXp": 1826,
    "chrXq": 2943,
    "chrYp": 5658,
    "chrYq": 6388,
}
