from teloBoundaryHelpers import *
from constants import expectedTeloCompositionQ, expectedTeloCompositionP, areaDiffsThreshold, teloNPTeloCompositionGStrand, teloNPTeloCompositionCStrand
import numpy as np
from Bio import SeqIO
import re


# The following function takes in a sequence, and returns the index of the telomere boundary.
def getTeloBoundary(seq, isGStrand = None, compositionGStrand=[], compositionCStrand = [], teloWindow=100, windowStep=6, changeThreshold=-20, plateauDetectionThreshold=-50, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=500, showGraphs=False, returnLastDiscontinuity=False, secondarySearch = False):
    """
    This function takes in a sequence, and returns the index of the telomere boundary.

    :param seq: The sequence to be analyzed
    :param isGStrand: True if the sequence is the G strand (has TTAGGG telomeres), False if it is the C strand (has CCCTAA telomeres). 
           Is None by default, and will be determined using the composition lists (count of G vs C repeats) if not specified.
    :param compositionGStrand: A list of lists, where each list contains a nucleotide pattern representing the expected telomere pattern on the G Strand.
    :param compositionCStrand: A list of lists, where each list contains a nucleotide pattern representing the expected telomere pattern on the C Strand.
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
    validate_seq_teloWindow(seq, teloWindow)

    if len(compositionGStrand) == 0:
        compositionGStrand = expectedTeloCompositionQ
    if len(compositionCStrand) == 0:
        compositionCStrand = expectedTeloCompositionP

    # calculate telomere strand type
    if isGStrand == None:
        isGStrand = getIsGStrandFromSeq(seq, compositionGStrand[targetPatternIndex], compositionCStrand[targetPatternIndex])
        if isGStrand < 0 or not isinstance(isGStrand, bool) and not isinstance(isGStrand, np.bool_):
            print("Could not determine telomere strand type, returning -1")
            return -1
        
    composition = []
    if isGStrand:
        composition = compositionGStrand
    else:
        composition = compositionCStrand


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
                # patternCount = the number of times we see the pattern
                # len(ntPattern) = the length of the pattern, so GGG is 3
                # teloLen = the length of the telomere window
                # patternComposition = the expected composition of the pattern, so 3/6 for GGG
                # Here we get an offset score based on the number of nucleotides we expect to see, given the pattern composition
                rawOffsetValue = (
                    (patternCount * len(ntPattern)) / teloLen - patternComposition)
                if patternComposition != 0:
                    # Here we convert the raw offset score to a percentage
                    percentOffsetValue = (
                        rawOffsetValue / patternComposition) * 100
                    currentOffsets.append(percentOffsetValue)
                else:
                    # Incase the expected composition is 0, we just return the raw offset score
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
        indexAtThreshold = (next((y for y in range(len(areaList) - 2, 0, -1) if (
            areaList[y] > changeThreshold and (0 > (areaList[y + 1] - areaList[y])))), indexAtThreshold))
        if indexAtThreshold != -1:
            # The min threshold was reached, and the slope was negative, so we can look for the max threshold ahead of it
            indexAtThreshold = (next((y for y in range(indexAtThreshold, len(areaList) - 2) if (
                areaList[y] < plateauDetectionThreshold and (0 > (areaList[y + 1] - areaList[y])))), indexAtThreshold))
        else:
            # Didn't find a point above the changeThreshold, so we just scan for the first point past the maxThreshold
            indexAtThreshold = (next((y for y in range(len(areaList) - 2) if (
                areaList[y] < plateauDetectionThreshold and (0 > (areaList[y + 1] - areaList[y])))), indexAtThreshold))
    else:
        indexAtThreshold = (next((y for y in range(len(areaList) - 2) if (
            areaList[y] < plateauDetectionThreshold and (0 > (areaList[y + 1] - areaList[y])))), indexAtThreshold))

    if indexAtThreshold == -1:
        print("No telo boundary found, returning -1")
        if showGraphs:
            graphLine(
                areaList, composition[targetPatternIndex][0] + " Area", windowStep)
            makeOffsetPlot(ntOffsets, composition,
                           offsetIndexToBPConstant=windowStep)
        return -1

    # Look through areaDiffs to find point where areaDiffs plateau
    for x in range(indexAtThreshold, len(areaDiffs)-1):
        # if it plateaus, or in the rare case that the diff jumps over the threshold, we have found the boundary point
        if abs(areaDiffs[x]) < areaDiffsThreshold or areaDiffs[x] < areaDiffsThreshold and areaDiffs[x+1] > areaDiffsThreshold :
            boundaryPoint = x * windowStep
            break
    if boundaryPoint == -1:
        # This means we have reached the end of the telomere
        # but we didn't fine the point at which the telomere offset stopped changing.
        print("Warning: Sequence was not long enough to find a telomere boundary, returning end of sequence as boundary point")
        boundaryPoint = len(areaDiffs) * windowStep
        
    if secondarySearch == True:
        # Lower being towards the telomere, upper being towards the centromere
        telomereOffset = 500
        subTelomereOffset = 1000
        telomereOffsetRE = 30
        subTelomereOffsetRE = 30
        ntPattern = ntPatternEntry[0]

        if isGStrand == True:
            # scanSeq = seq[boundaryPoint-telomereOffset:boundaryPoint+subTelomereOffset]
            lowerIndex = (len(seq) - (boundaryPoint+subTelomereOffset))
            upperIndex = len(seq) - (boundaryPoint-telomereOffset)            
            if lowerIndex < 0:
                lowerIndex = 0
            if upperIndex > len(seq):
                upperIndex = len(seq)
            
            scanSeq = seq[lowerIndex:upperIndex]
            if len(scanSeq) < 100:
                print("Warning: Sequence was not long enough to perform secondary search, returning original boundary point")
            else:
                secBoundary = getTeloBoundary(scanSeq, isGStrand, composition, teloWindow=90, windowStep=6, changeThreshold=changeThreshold, plateauDetectionThreshold=-60, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=100, showGraphs=showGraphs, returnLastDiscontinuity=returnLastDiscontinuity, secondarySearch=False)

                tempBoundary = boundaryPoint + secBoundary - telomereOffset
                if upperIndex == len(seq):
                    tempBoundary = secBoundary

                # ***( upper and lower here might be wrong, check this)
                scanSeq = seq[(len(seq) - (tempBoundary+subTelomereOffsetRE)):len(seq) - (tempBoundary-telomereOffsetRE)]
                scan_pattern = "("+ntPattern+")" + "("+ntPattern+")"
                match = re.search(str(scan_pattern), str(scanSeq))
                if match:
                    secBoundary = len(scanSeq) -match.span()[0]
                    boundaryPoint = tempBoundary + secBoundary - telomereOffsetRE
                else:
                    boundaryPoint = tempBoundary
                    print("Secondary search failed to find a match, returning original boundary point")
        else:
            # print("Performing secondary search")
            lowerIndex = boundaryPoint-telomereOffset
            upperIndex = boundaryPoint+subTelomereOffset
            if lowerIndex < 0:
                lowerIndex = 0
            if upperIndex > len(seq):
                upperIndex = len(seq)
            scanSeq = seq[lowerIndex:upperIndex]
            if len(scanSeq) < 100:
                print("Warning: Sequence was not long enough to perform secondary search, returning original boundary point")
            else:
                secBoundary = getTeloBoundary(scanSeq, isGStrand, composition, teloWindow=90, windowStep=6, changeThreshold=changeThreshold, plateauDetectionThreshold=-60, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=100, showGraphs=showGraphs, returnLastDiscontinuity=returnLastDiscontinuity, secondarySearch=False)

                tempBoundary = boundaryPoint + secBoundary - telomereOffset
                if lowerIndex == 0:
                    tempBoundary = secBoundary
                scanSeq = seq[tempBoundary-telomereOffsetRE:tempBoundary+subTelomereOffsetRE]
                
                if patternComposition != 1:
                    print("Warning: Secondary search is not fully compatible with telomere compositions less than 1. Please provide a telomere pattern that covers 6/6 of the expected telomere nucleotides, like 'TTAGGG' or 'GGG...'.")
                else:
                    
                    scan_pattern = "("+ntPattern+")" + "("+ntPattern+")"
                    matches = [match for match in re.finditer(scan_pattern, str(scanSeq))]
                    if matches:
                        teloEnd = matches[-1].end()
                        boundaryPoint = tempBoundary + teloEnd - telomereOffsetRE
                    else:
                        boundaryPoint = tempBoundary
                        print("Secondary search failed to find a match, returning original boundary point")

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

def trimTeloReferenceGenome(filename, outputFilename, subSec = None, compositionCStrandIn=expectedTeloCompositionP, compositionGStrandIn=expectedTeloCompositionQ, teloWindowIn=100, windowStepIn=6, plateauDetectionThresholdIn=-60, changeThresholdIn=-20, targetPatternIndexIn=-1, nucleotideGraphAreaWindowSizeIn=500, showGraphsIn=False, returnLastDiscontinuityIn=False, secondarySearchIn = False):
    trimmed_sequences = []

    # Trims the records and saves them
    for record in SeqIO.parse(filename, "fasta"):
        startTeloLength = getTeloBoundary(record.seq[:500000], isGStrand=False,  compositionGStrand = compositionGStrandIn, compositionCStrand=compositionCStrandIn, teloWindow=teloWindowIn, windowStep=windowStepIn, plateauDetectionThreshold=plateauDetectionThresholdIn, changeThreshold=changeThresholdIn,
                                          targetPatternIndex=targetPatternIndexIn, nucleotideGraphAreaWindowSize=nucleotideGraphAreaWindowSizeIn, showGraphs=showGraphsIn, returnLastDiscontinuity=returnLastDiscontinuityIn, secondarySearch = secondarySearchIn)
        endTeloLength = getTeloBoundary(record.seq[-500000:], isGStrand=True,  compositionGStrand = compositionGStrandIn, compositionCStrand=compositionCStrandIn, teloWindow=teloWindowIn, windowStep=windowStepIn, plateauDetectionThreshold=plateauDetectionThresholdIn, changeThreshold=changeThresholdIn,
                                        targetPatternIndex=targetPatternIndexIn, nucleotideGraphAreaWindowSize=nucleotideGraphAreaWindowSizeIn, showGraphs=showGraphsIn, returnLastDiscontinuity=returnLastDiscontinuityIn, secondarySearch = secondarySearchIn)
        if subSec == None:
            trimmed_sequences.append(record[startTeloLength:-endTeloLength])
        else:
            # We will save 2 subsequences, one for the q arm and one for the p arm
            cStrand = record[startTeloLength:startTeloLength+subSec]
            gStrand = record[-endTeloLength-subSec:-endTeloLength]
            # change the id of the sequence
            cStrand.id = cStrand.id + "_C"
            gStrand.id = gStrand.id + "_G"
            trimmed_sequences.append(cStrand)
            trimmed_sequences.append(gStrand)

    # Write the trimmed sequences to the output file
    SeqIO.write(trimmed_sequences, outputFilename, "fasta")

def getIsGStrandFromSeq(seq, GStrandPatternIn, CStrandPatternIn, searchStrandRepeats = 4, minTeloCountDiff = 1, fusedReadTeloRepeatThreshold = 20, maxTelomereGap = 150):
    # We will look at the beginning of the seq and count for C strands, then look at the end and count for G strands
    # the compare the counts to see which is greater and return the result
    
    CStrandPattern = CStrandPatternIn[0]
    GStrandPattern = GStrandPatternIn[0]

    boundaryReg = ".{0,6}"

    CStrandPattern =  ("("+CStrandPattern+")"+boundaryReg) * searchStrandRepeats
    GStrandPattern =  ("("+GStrandPattern+")"+boundaryReg) * searchStrandRepeats
    CStrandPattern = CStrandPattern[:-len(boundaryReg)]
    GStrandPattern = GStrandPattern[:-len(boundaryReg)]

    allCStrands = [match for match in re.finditer(CStrandPattern, str(seq.upper()))]
    # get start index of each match. For C strand start at the beginning and count forward
    cStrandCount = 0
    cStrandMatchLengths = 0
    lastMatch = 0
    for match in allCStrands:
        if lastMatch == 0:
            lastMatch = match.start()

        if match.start() - lastMatch <= maxTelomereGap:
            cStrandCount += 1
            cStrandMatchLengths += match.end() - match.start()
        lastMatch = match.start()
        
    allGStrands = [match for match in re.finditer(GStrandPattern, str(seq.upper()))]

    gStrandCount = 0
    gStrandMatchLengths = 0
    lastMatch = 0
    for match in allGStrands[::-1]:
        if lastMatch == 0:
            lastMatch = match.start()
        if lastMatch - match.start() <= maxTelomereGap:
            gStrandCount += 1
            gStrandMatchLengths += match.end() - match.start()
        lastMatch = match.start()

    if max(cStrandCount, gStrandCount) <= minTeloCountDiff and abs(cStrandCount - gStrandCount) <= minTeloCountDiff:

        # The intention here is to give a last chance to classify the read, as we want to keep as many short reads as possible
        if cStrandMatchLengths > gStrandMatchLengths:
            return False
        elif gStrandMatchLengths > cStrandMatchLengths:
            return True
        print("Warning: could not determine telomere strand type from sequence, returning -20")
        return -20
    
    if min(cStrandCount, gStrandCount) > 100 or abs(cStrandCount - gStrandCount) <= 0.6 * min(cStrandCount, gStrandCount) or min(cStrandCount, gStrandCount) >= fusedReadTeloRepeatThreshold:
        print("Warning: fused strand likely, returning -10")
        return -10

    if cStrandCount > gStrandCount:
        return False
    else:
        return True


def isGStrand(chrArm,strand):
    if strand == "+" and chrArm == "q":
        return True
    elif strand == "-" and chrArm == "p":
        return True
    elif strand == "+" and chrArm == "p":
        return False
    elif strand == "-" and chrArm == "q":
        return False
    else:
        print("error, could not identify strand")
        return False
    


def getTeloNPBoundary(seq, isGStrand=None, compositionCStrandIn=teloNPTeloCompositionCStrand, compositionGStrandIn=teloNPTeloCompositionGStrand, teloWindow=100, windowStep=6, changeThreshold=-20, plateauDetectionThreshold=-60, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=750, showGraphs=False, returnLastDiscontinuity=True, secondarySearch = True):
    return getTeloBoundary(seq, isGStrand, compositionCStrand=compositionCStrandIn, compositionGStrand=compositionGStrandIn, teloWindow=teloWindow, windowStep=windowStep, changeThreshold=changeThreshold, plateauDetectionThreshold=plateauDetectionThreshold, targetPatternIndex=targetPatternIndex, nucleotideGraphAreaWindowSize=nucleotideGraphAreaWindowSize, showGraphs=showGraphs, returnLastDiscontinuity=returnLastDiscontinuity, secondarySearch=secondarySearch)
