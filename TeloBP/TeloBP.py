from teloBoundaryHelpers import *
from constants import expectedTeloCompositionQ, expectedTeloCompositionP, areaDiffsThreshold
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
# import matplotlib.pyplot as plt


# This value must be large enough that it can "look into" the area well past the telomere boundary, allowing it to
# make sure that the discontinuity in the telomere pattern is sustained, and that the first change in the pattern
# can be identified. If it is too small, then the first change in the area graph may not be idenfied, as the
# offset values will be too noisy.
# nucleotideGraphAreaWindowSize = 500
# nucleotideGraphAreaWindowSize = 10


# I don't know how I can properly explain the areaThreshold values, without just showing you the graphs
def getTeloBoundary(seq, isGStrand, composition=[], teloWindow=100, windowStep=6, changeThreshold=-20, plateauDetectionThreshold=-50, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=500, showGraphs=False, returnLastDiscontinuity=False):
    """
    This function takes in a sequence, and returns the index of the telomere boundary.

    :param seq: The sequence to be analyzed
    :param isGStrand: True if the sequence is the G strand (has TTAGGG telomeres), False if it is the C strand (has CCCTAA telomeres). 
    :param composition: A list of lists, where each list contains a nucleotide pattern, and the expected composition of that pattern 
           in the telomere.
    :param teloWindow: The size of the window to be used when calculating the offset of the nucleotide composition from the expected 
           telomere composition.
    :param windowStep: The step size to be used when moving through the sequence in 'windows' of size teloWindow.
    :param changeThreshold: The sequence difference threshold where we can assume that we are approaching the telomere boundary. 
           This value is only used if returnLastDiscontinuity is true. 
    :param plateauDetectionThreshold: The sequence difference threshold at which the algorithm starts scanning for the point where the 
           slope of discontinuity plateaus. In this context, discontinuity means when the sequence no longer has a pattern similar to a telomere. 
    :param targetPatternIndex: The index of the pattern in the composition list that we want to use to calculate the telomere boundary. 
           This is primarily for testing purposes, when we want to visualize the offsets of multiple patterns, but only use 1 for finding the boundary.
    :param nucleotideGraphAreaWindowSize: The size of the window to be used when calculating the area under the curve of the offsets. Generally speaking, 
           a smaller value gives greater precision, but less accuracy if the inputted sequence is noisy. A larger value gives greater accuracy, but less precision.
           More specifically, this value must be large enough that it can "look into" the area well past the telomere boundary, making sure that the 
           discontinuity in the telomere pattern is sustained.
    :param showGraphs: Boolean value, if true, the graphs will be shown.
    :param returnLastDiscontinuity: Boolean value, if true, the algorithm will return the last possible point of discontinuity, rather than the first.
           For sequences which are very noisy, this may be necessary, as the first point of discontinuity may be a false positive.
    """

    boundaryPoint = -1
    ntOffsets = []
    graphAreaWindowSize = int(nucleotideGraphAreaWindowSize / windowStep)

    if len(composition) == 0:
        # print("Warning: no composition list provided, using default telomere compositions")
        if isGStrand == True:
            composition = expectedTeloCompositionQ
        else:
            composition = expectedTeloCompositionP

    validate_parameters(seq, isGStrand, composition, teloWindow, windowStep, plateauDetectionThreshold,
                        changeThreshold, targetPatternIndex, nucleotideGraphAreaWindowSize, showGraphs)

    # Move through the sequence in windows of size teloWindow, and step size windowStep,
    # and calculate the offset of the nucleotide composition from the expected telomere composition
    for i in range(0, len(seq) - teloWindow, windowStep):
        teloSeq = ""
        if isGStrand == True:
            teloSeq = seq[(len(seq) - i) - teloWindow:len(seq) - i]
        else:
            teloSeq = seq[i:i + teloWindow]
        teloLen = len(teloSeq)
        teloSeqUpper = str(teloSeq.upper())

        currentOffsets = []
        # Calculate the offset of the nucleotide composition from the expected telomere composition
        for ntPatternEntry in composition:
            ntPattern = ntPatternEntry[0]
            patternComposition = ntPatternEntry[1]
            patternCount = len(re.findall(ntPattern, teloSeqUpper))
            if is_regex_pattern(ntPattern) == True:
                if len(ntPatternEntry) != 3:
                    print(
                        "Error: a target length must be specified as a third list item if using a regex pattern. Example: ['GGG|AAA', 3/6, 3], where the third item is the target length.")
                    return -1
                regexTargetLength = ntPatternEntry[2]
                rawOffsetValue = (
                    (patternCount * regexTargetLength) / teloLen - patternComposition)
                if patternComposition != 0:
                    percentOffsetValue = (
                        rawOffsetValue / patternComposition) * 100
                    currentOffsets.append(percentOffsetValue)
                else:
                    percentOffsetValue = (rawOffsetValue) * 100
                    currentOffsets.append(percentOffsetValue)
            else:
                rawOffsetValue = ((patternCount * len(ntPattern)) /
                                  teloLen - patternComposition)
                if patternComposition != 0:
                    percentOffsetValue = (
                        rawOffsetValue / patternComposition) * 100
                    currentOffsets.append(percentOffsetValue)
                else:
                    percentOffsetValue = (rawOffsetValue) * 100
                    currentOffsets.append(percentOffsetValue)

        ntOffsets.append(currentOffsets)

    areaList = getGraphArea(ntOffsets, targetPatternIndex, graphAreaWindowSize)
    areaDiffs = np.diff(areaList)
    indexAtThreshold = -1

    # If returnLastDiscontinuity is true, we will scan for the
    # last point where we are above the changeThreshold, then look ahead for the first point where we
    # go below the plateauDetectionThreshold. This is because the area under the curve is not always monotonically
    # decreasing.

    if returnLastDiscontinuity:
        # Here, we grab the last point where the area is below the changeThreshold, and the slope is negative
        # *** 0<(areaList[y]-areaList[y+1] is a weird way of checking if the slope is negative, but it works. might change later
        indexAtThreshold = (next((y for y in range(len(areaList) - 2, 0, -1) if (
            areaList[y] > changeThreshold and (0 < (areaList[y] - areaList[y + 1])))), indexAtThreshold))
        if indexAtThreshold != -1:
            # The min threshold was reached, and the slope was negative, so we can look for the max threshold ahead of it
            indexAtThreshold = (next((y for y in range(indexAtThreshold, len(areaList) - 2) if (
                areaList[y] < plateauDetectionThreshold and (0 < (areaList[y] - areaList[y + 1])))), indexAtThreshold))
        else:
            # Didn't find a point above the changeThreshold, so we just scan for the first point past the maxThreshold
            indexAtThreshold = (next((y for y in range(len(areaList) - 2) if (
                areaList[y] < plateauDetectionThreshold and (0 < (areaList[y] - areaList[y + 1])))), indexAtThreshold))
    else:
        indexAtThreshold = (next((y for y in range(len(areaList) - 2) if (
            areaList[y] < plateauDetectionThreshold and (0 < (areaList[y] - areaList[y + 1])))), indexAtThreshold))

    if indexAtThreshold == -1:
        print("No telo boundary found, returning -1")
        if showGraphs:
            graphLine(
                areaList, composition[targetPatternIndex][0] + " Area", windowStep)
            makeOffsetPlot(ntOffsets, composition,
                           offsetIndexToBPConstant=windowStep)
        return -1

    # Look through areaDiffs to find point where areaDiffs plateau
    for x in range(indexAtThreshold, len(areaDiffs)):
        if abs(areaDiffs[x]) < areaDiffsThreshold:
            boundaryPoint = x * windowStep
            break
    if boundaryPoint == -1:
        # This means we have reached the end of the telomere
        # but we didn't fine the point at which the telomere offset stopped changing.
        print("Warning: Sequence was not long enough to find a telomere boundary, returning end of sequence as boundary point")
        boundaryPoint = len(areaDiffs) * windowStep

    if showGraphs:
        graphLine(areaList, composition[targetPatternIndex]
                  [0] + " Area", windowStep, boundaryPoint=boundaryPoint)
        makeOffsetPlot(ntOffsets, composition, windowStep)

    return boundaryPoint


# The following function takes in the location of a reference genome, and outputs
# a genome file with the telomeres removed. The teloBP parameters have been set
# to the recommended values mentioned in the Readme, but can be modified.
# NOTE: This algorithm assumes that each chromosome is its own read, and not split into
# q and p arms.

def trimTeloReferenceGenome(filename, outputFilename, compositionIn=[], teloWindowIn=100, windowStepIn=6, plateauDetectionThresholdIn=-60, changeThresholdIn=-20, targetPatternIndexIn=-1, nucleotideGraphAreaWindowSizeIn=500, showGraphsIn=False, returnLastDiscontinuityIn=False):
    trimmed_sequences = []

    # Trims the records and saves them
    for record in SeqIO.parse(filename, "fasta"):
        startTeloLength = getTeloBoundary(record.seq[:500000], isGStrand=False, composition=compositionIn, teloWindow=teloWindowIn, windowStep=windowStepIn, plateauDetectionThreshold=plateauDetectionThresholdIn, changeThreshold=changeThresholdIn,
                                          targetPatternIndex=targetPatternIndexIn, nucleotideGraphAreaWindowSize=nucleotideGraphAreaWindowSizeIn, showGraphs=showGraphsIn, returnLastDiscontinuity=returnLastDiscontinuityIn)
        endTeloLength = getTeloBoundary(record.seq[-500000:], isGStrand=True, composition=compositionIn, teloWindow=teloWindowIn, windowStep=windowStepIn, plateauDetectionThreshold=plateauDetectionThresholdIn, changeThreshold=changeThresholdIn,
                                        targetPatternIndex=targetPatternIndexIn, nucleotideGraphAreaWindowSize=nucleotideGraphAreaWindowSizeIn, showGraphs=showGraphsIn, returnLastDiscontinuity=returnLastDiscontinuityIn)
        trimmed_sequences.append(record[startTeloLength:-endTeloLength])

    # Write the trimmed sequences to the output file
    SeqIO.write(trimmed_sequences, outputFilename, "fasta")
