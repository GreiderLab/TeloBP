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
#     ["T", expectedTeloCompositionQ[1][1]],
#     ["A", expectedTeloCompositionQ[0][1]],
#     ["G", expectedTeloCompositionQ[3][1]],
#     ["C", expectedTeloCompositionQ[2][1]],
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

# Recommended composition for nanopore reads called with Guppy
teloNPTeloCompositionGStrand = [["(?!GGG)[ATC]{2}.GGG|(?!AAA)[AT]{3}AAA|TTAGG.", 6/6, 6]]

teloNPTeloCompositionCStrand = [["CTTCTT|CCTGG|CCC(?!CCC)[ATCG]{3}", 6/6, 6]]

# Recommended composition for R10 nanopore reads: the canonical repeat plus every 1 bp substitution of it
teloNPR10TeloCompositionGStrand = [["TTAGGG|[ACG]TAGGG|T[ACG]AGGG|TT[CGT]GGG|TTA[ACT]GG|TTAG[ACT]G|TTAGG[ACT]", 6/6, 6]]

teloNPR10TeloCompositionCStrand = [["CCCTAA|[AGT]CCTAA|C[AGT]CTAA|CC[AGT]TAA|CCC[ACG]AA|CCCT[CGT]A|CCCTA[CGT]", 6/6, 6]]

# Recommended composition for R9 nanopore reads called with Guppy: the R10 patterns plus Guppy's telomere miscalls.
# The Guppy patterns alone miss some variant repeats common near the subtelomere, such as TGAGGG.
teloNPR9TeloCompositionGStrand = [[teloNPR10TeloCompositionGStrand[0][0] + "|" + teloNPTeloCompositionGStrand[0][0], 6/6, 6]]

teloNPR9TeloCompositionCStrand = [[teloNPR10TeloCompositionCStrand[0][0] + "|" + teloNPTeloCompositionCStrand[0][0], 6/6, 6]]

# getTeloNPBoundary settings for each nanopore chemistry.
# R10: for R10.4.1 reads (Dorado). The first discontinuity is returned, stepping over gaps only when
#      telomere sequence resumes within maxDiscontinuityGap bp, so a distant ITS is not counted as telomere.
# R9: for R9.4.1 reads basecalled with Guppy. The patterns accept every R10 repeat and Guppy's telomere
#     miscalls, with the same bounded gap search as R10.
# R9_archive: the original TeloNP settings (v1), kept unchanged so earlier results can be reproduced.
#     It returns the last discontinuity in the read, so an ITS far from the telomere can be counted as telomere.
teloNPPresets = {
    "R9": {
        "compositionGStrand": teloNPR9TeloCompositionGStrand,
        "compositionCStrand": teloNPR9TeloCompositionCStrand,
        "nucleotideGraphAreaWindowSize": 750,
        "returnLastDiscontinuity": False,
        "maxDiscontinuityGap": 2000,
    },
    "R9_archive": {
        "compositionGStrand": teloNPTeloCompositionGStrand,
        "compositionCStrand": teloNPTeloCompositionCStrand,
        "nucleotideGraphAreaWindowSize": 750,
        "returnLastDiscontinuity": True,
        "maxDiscontinuityGap": 0,
    },
    "R10": {
        "compositionGStrand": teloNPR10TeloCompositionGStrand,
        "compositionCStrand": teloNPR10TeloCompositionCStrand,
        "nucleotideGraphAreaWindowSize": 750,
        "returnLastDiscontinuity": False,
        "maxDiscontinuityGap": 2000,
    },
}


# When the area under the curve of the offsets is calculated, this
# constant is used to mark the point where the slope of the curve
# starts to plateau. Basically, once the difference between two
# values in the area under the curve is less than this, the slope
# is considered to be plateauing.
areaDiffsThreshold = 0.2



errorReturns = {"init": -1, "fusedRead": -10,"strandType": -20, "seqNotFound": -1000}


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
