from teloBoundaryHelpers import *
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
# import matplotlib.pyplot as plt
import re


# This new version of the program will take in an expected telomere composition dictionary, and will 
# also take in a threshold at which it will consider its current window no longer telomere, and return the telomere
# boundary. 
# "TTAGGG"


# by default the length of the pattern is the length value, but then it 
# needs to be specified for the regex patterns.
expectedTeloCompositionQ = [
    ["T", 2/6],
    ["A", 1/6],
    ["G", 3/6],
    ["C", 0/6],
    ["GGG", 3/6],
    ["TT", 2/6],
    ["GGG|AAA", 3/6, 3],
]


expectedTeloCompositionP = [
    ["T", expectedTeloCompositionQ[0][1]],
    ["A", expectedTeloCompositionQ[1][1]],
    ["G", expectedTeloCompositionQ[2][1]],
    ["C", expectedTeloCompositionQ[3][1]],
    ["CCC", expectedTeloCompositionQ[4][1]],
    ["AA", expectedTeloCompositionQ[5][1]],
    ["CCCTAA|CTTCTT|CCCTGG|CCTGG", 6/6, 6],
    ["CCC|TTT", expectedTeloCompositionQ[6][1],3]
]


areaDiffsThreshold = 0.1
# This value must be large enough that it can "look into" the area well past the telomere boundary, allowing it to 
# make sure that the discontinuity in the telomere pattern is sustained, and that the first change in the pattern
# can be identified. If it is too small, then the first change in the area graph may not be idenfied, as the 
# offset values will be too noisy. 
# nucleotideGraphAreaWindowSize = 500
# nucleotideGraphAreaWindowSize = 10

def is_regex_pattern(input_string):
    return bool(re.search(r'[^a-zA-Z]', input_string))

# I don't know how I can properly explain the areaThreshold values, without just showing you the graphs
def getTeloBoundary(seq, isGStrand, composition = [], teloWindow = 100, windowStep = 6, maxAreaThreshold = -15, minAreaThreshold = -5, targetPatternIndex=-1, nucleotideGraphAreaWindowSize = 500, showGraphs = False, returnLastDiscontinuity = False):
    
    """
    This function takes in a sequence, and returns the index of the telomere boundary.

    :param seq: The sequence to be analyzed
    :param isGStrand: True if the sequence is the G strand (has TTAGGG telomeres), False if it is the C strand (has CCCTAA telomeres). 
    :param composition: A list of lists, where each list contains a nucleotide pattern, and the expected composition of that pattern in the telomere.
    :param teloWindow: The size of the window to be used when calculating the offset of the nucleotide composition from the expected telomere composition.
    :param windowStep: The step size to be used when moving through the sequence in 'windows' of size teloWindow.
    :param maxAreaThreshold: The maximum value that the area under the curve can reach before the program will stop looking for the telomere boundary.
    """
    
    # print(seq)
    boundaryPoint = -1
    ntOffsets = []
    graphAreaWindowSize = int(nucleotideGraphAreaWindowSize/windowStep)
    # print("graphAreaWindowSize: " + str(graphAreaWindowSize))

    if len(composition) == 0:
        print("Warning: no composition list provided, using default telomere compositions")
        if isGStrand == True:
            composition = expectedTeloCompositionQ
        else:
            # print("c strand")
            composition = expectedTeloCompositionP

    validate_parameters(seq, isGStrand, composition, teloWindow, windowStep, maxAreaThreshold, minAreaThreshold, targetPatternIndex, nucleotideGraphAreaWindowSize, showGraphs)
    
    # Move through the sequence in windows of size teloWindow, and step size windowStep, 
    # and calculate the offset of the nucleotide composition from the expected telomere composition
    for i in range(0, len(seq)-teloWindow, windowStep):
        teloSeq = ""
        if isGStrand == True:
            teloSeq = seq[(len(seq)-i)-teloWindow:len(seq)-i]
        else:
            teloSeq = seq[i:i+teloWindow]
        teloLen = len(teloSeq)
        teloSeqUpper = str(teloSeq.upper())
        seq_object = Seq(teloSeq.upper())

        currentOffsets = []
        # Calculate the offset of the nucleotide composition from the expected telomere composition
        for ntPatternEntry in composition:
            ntPattern = ntPatternEntry[0]
            patternComposition = ntPatternEntry[1]
            patternCount = len(re.findall(ntPattern, teloSeqUpper))
            if is_regex_pattern(ntPattern) == True:
                if len(ntPatternEntry) != 3:
                    print("Error: a target length must be specified as a third list item if using a regex pattern. Example: ['GGG|AAA', 3/6, 3], where the third item is the target length.")
                    return -1
                regexTargetLength = ntPatternEntry[2]
                rawOffsetValue = ((patternCount*regexTargetLength)/teloLen - patternComposition)
                if patternComposition !=0:
                    percentOffsetValue = (rawOffsetValue/patternComposition) * 100
                    currentOffsets.append(percentOffsetValue)
                else:
                    percentOffsetValue = (rawOffsetValue) * 100
                    currentOffsets.append(percentOffsetValue)
            else:
                rawOffsetValue = ((patternCount*len(ntPattern))/teloLen - patternComposition)
                if patternComposition !=0:
                    percentOffsetValue = (rawOffsetValue/patternComposition) * 100
                    currentOffsets.append(percentOffsetValue)
                else:
                    percentOffsetValue = (rawOffsetValue) * 100
                    currentOffsets.append(percentOffsetValue)
            
            # currentOffsets.append((counts[ntPattern])/teloLen - composition[ntPattern])
            # This one gives us a percentage offset, so 0 is no offset, -1 is 100% below expected, 1 is 100% above expected

            # ******* come back to this later and figure out a new offset scoring calculation


            # currentOffsets.append(((counts[ntPattern])/teloLen)/composition[ntPattern] - 1)
        ntOffsets.append(currentOffsets)
        

    # print(ntOffsets)
    areaList = getGraphArea(ntOffsets, targetPatternIndex, graphAreaWindowSize)
    # print(areaList)
    areaDiffs = np.diff(areaList)
    indexAtThreshold = -1

    newAreaDiffs = areaList.copy()
    # print(newAreaDiffs)

    # try:
    # before I would scan for the point where we first go over the threshold, but now I am scanning for the 
    # last point where we are bellow the minAreaThreshold, and then I will scan for the first point where we 
    # go above the maxAreaThreshold. This is because the area under the curve is not always monotonically
    # decreasing.

    if returnLastDiscontinuity:
        # Here, we grab the last point where the area is below the minAreaThreshold, and the slope is negative
        # *** 0<(areaList[y]-areaList[y+1] is a weird way of checking if the slope is negative, but it works. might change later
        indexAtThreshold = (next((y for y in range(len(areaList)-2, 0, -1) if (areaList[y] > minAreaThreshold and (0<(areaList[y]-areaList[y+1]) ))), indexAtThreshold))
        if indexAtThreshold != -1 :
            # The min threshold was reached, and the slope was negative, so we can look for the max threshold ahead of it
            indexAtThreshold = (next((y for y in range(indexAtThreshold,len(areaList)-2) if (areaList[y] < maxAreaThreshold and (0<(areaList[y]-areaList[y+1]) ))), indexAtThreshold))
        else:
            # Didn't find a point above the minAreaThreshold, so we just scan for the first point past the maxThreshold
            indexAtThreshold = (next((y for y in range(len(areaList)-2) if (areaList[y] < maxAreaThreshold and (0<(areaList[y]-areaList[y+1]) ))), indexAtThreshold))
    else:
        indexAtThreshold = (next((y for y in range(len(areaList)-2) if (areaList[y] < maxAreaThreshold and (0<(areaList[y]-areaList[y+1]) ))), indexAtThreshold))


    if indexAtThreshold == -1:
        print("No telo boundary found on q end")
        if showGraphs:
            # print("First show graph called")
            graphLine(areaList, composition[targetPatternIndex][0]+" Area", windowStep)
            makeOffsetPlot(ntOffsets,composition, offsetIndexToBPConstant= windowStep)
        return -1

    # print("boundary index start: ", (indexAtThreshold * windowStep))
    # Look through areaDiffs to find first index with difference close to zero
    for x in range(indexAtThreshold, len(areaDiffs)):
        if abs(areaDiffs[x]) < areaDiffsThreshold:
            # print("passed the check")
            # print(areaDiffs[x-5:x+5])
            # print("indexAtThreshold: ", indexAtThreshold)
            # print("boundaryPoint: ", x)
            # print("graphArea window: ", graphAreaWindowSize)
            boundaryPoint = x * windowStep
            break
    if boundaryPoint == -1:
        # This means we have reached the threshold value indicating that we have reached the end of the telomere
        # but we didn't fine the point at which the telomere offset stopped changing. 
        print("Didn't find the point at which the telomere offset stopped changing.")
        boundaryPoint = len(areaDiffs) * windowStep

    if showGraphs:
        # print("second show graph called")
        graphLine(areaList, composition[targetPatternIndex][0]+" Area", windowStep, boundaryPoint = boundaryPoint)
        makeOffsetPlot(ntOffsets,composition, windowStep)

    # print("q returning ", boundaryPoint)
    return boundaryPoint
        