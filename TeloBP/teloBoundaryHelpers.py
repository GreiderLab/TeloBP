import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

from TeloBP.constants import manualLabelsCHM13Positions, manualLabelsCHM13, areaDiffsThreshold
import warnings


def descriptionToChr(description):
    # print(description)
    chrNum = description.split(" ")[-1]
    # if chrnum is a number
    if not chrNum.isdigit():
        if chrNum == "X" or chrNum == "Y":
            return "chr_" + chrNum
        return chrNum
    if len(chrNum) == 1:
        output = "chr_0" + chrNum
        return output
    else:
        output = "chr_" + chrNum
        return output


def descriptionToChrName(description):
    return description.split(" ")[0]


def is_regex_pattern(input_string):
    return bool(re.search(r'[^a-zA-Z]', input_string))


def write_bed_file(file_path, bed_data):
    with open(file_path, 'w') as bed_file:
        for entry in bed_data:
            bed_file.write('\t'.join(str(e) for e in entry) + '\n')


def getGraphArea(offsets, targetColumn, windowSize):
    data = np.array(offsets)
    transposed_data = data.T
    areaList = []
    row = transposed_data[targetColumn, :]
    
    # print(row)
    for i in range(0, len(row) - windowSize, 1):
        area = row[i:i + windowSize].sum()
        areaList.append((area / windowSize))
        # areaList.append(area)
    return areaList


def getPlateauIndex(areaDiffs, startIndex):
    # Look through areaDiffs to find point where areaDiffs plateau. Returns -1 if they never do.
    for x in range(startIndex, len(areaDiffs)-1):
        # if it plateaus, or in the rare case that the diff jumps over the threshold, we have found the boundary point
        if abs(areaDiffs[x]) < areaDiffsThreshold or areaDiffs[x] < areaDiffsThreshold and areaDiffs[x+1] > areaDiffsThreshold :
            return x
    return -1


def stepOverShortGaps(areaList, areaDiffs, plateauIndex, maxGapSteps, changeThreshold, plateauDetectionThreshold):
    # If the area climbs back above changeThreshold within maxGapSteps of the boundary, the dip was a gap
    # inside the telomere (e.g. a stretch of bad basecalls), so we look for the next discontinuity instead.
    # Telomeric sequence further away than that, like an ITS, is left out of the telomere.
    while plateauIndex != -1:
        resumeIndex = next((y for y in range(plateauIndex, min(len(areaList), plateauIndex + maxGapSteps))
                            if areaList[y] > changeThreshold), -1)
        if resumeIndex == -1:
            return plateauIndex
        nextIndex = next((y for y in range(resumeIndex, len(areaList) - 2) if (
            areaList[y] < plateauDetectionThreshold and (0 > (areaList[y + 1] - areaList[y])))), -1)
        if nextIndex == -1:
            # The telomere continues to the end of the sequence
            return -1
        plateauIndex = getPlateauIndex(areaDiffs, nextIndex)
    return plateauIndex


def graphLine(rowIn, labelIn, windowStep, boundaryPoint=-1, pdfOut=None):

    row = rowIn
    # if boundaryPoint is less than half the x axis, only show boundaryPoint x2 of the data
    if boundaryPoint < 0.5 * len(rowIn):
        row = rowIn[:int((2 * boundaryPoint)/windowStep)]

    x = np.arange(len(row))
    x = x * windowStep
    fig, ax = plt.subplots()

    ax.plot(x, row, label=labelIn)

    # Set labels and title
    ax.set_xlabel('Distance from end of sequence')
    ax.set_ylabel('Nucleotide offset')
    ax.set_title('Nucleotide Offsets from Expected Telomere Composition')

    if boundaryPoint != -1:
        ax.axvline(boundaryPoint, color='red', label='Boundary Point')

    # Show legend
    ax.legend()

    if pdfOut is not None:
        pdfOut.savefig(fig)
    elif pdfOut is None:
        # Display the graph
        plt.show()


def makeOffsetPlot(offsets, compositions, offsetIndexToBPConstant):
    data = np.array(offsets)
    transposed_data = data.T

    x = np.arange(transposed_data.shape[1])

    # multiple each value in x by offsetIndexToBPConstant
    x = x * offsetIndexToBPConstant
    fig, ax = plt.subplots()

    labelDict = {}
    for patternI in range(len(compositions)):
        pattern = compositions[patternI][0]
        labelDict[patternI] = pattern

    # Plot each line
    for i in range(len(labelDict)):
        # This line reverses the order of the data. Used for figure in the paper
        # ax.plot(x, transposed_data[i, ::-1], label=f"{labelDict[i]} Offset")

        ax.plot(x, transposed_data[i, :], label=f"{labelDict[i]} Offset")

    # Set labels and title
    ax.set_xlabel('Distance from end of sequence (bp)')
    ax.set_ylabel('Nucleotide offset (%)')
    ax.set_title('Nucleotide Offsets from Expected Telomere Composition')
    ax.legend()
    plt.show()


def testTeloGenomePosition(chr, pos, testDict=manualLabelsCHM13Positions):
    # print(chr)
    if chr not in testDict:
        print("Error: chromosome not found in dictionary, returning 0")
        return 0
    expectedPos = testDict[chr]
    # print("Expected position: " + str(expectedPos))
    print(chr + " offset: " + str(pos - expectedPos) + "bp (obs - exp)")
    return (pos - expectedPos)


def testTeloLength(chr, length, testDict=manualLabelsCHM13):
    if chr not in testDict:
        print("Error: chromosome not found in dictionary, returning 0")
        return 0
    expectedLength = testDict[chr]
    # print(chr + " offset: " + str(length - expectedLength) + "bp (obs - exp)")
    return (length - expectedLength)

def recordBedData(bed_data, record, startTeloLength, endTeloLength):
    chrName = record.id
    pEndIndex = startTeloLength
    bed_data.append([chrName, pEndIndex-6, pEndIndex, chrName])
    
    qEndIndex = len(record.seq)-endTeloLength
    bed_data.append([chrName, qEndIndex, qEndIndex+6, chrName])

def write_bed_file(file_path, bed_data):
    with open(file_path, 'w') as bed_file:
        for entry in bed_data:
            bed_file.write('\t'.join(str(e) for e in entry) + '\n')


def validate_parameters(seq, isGStrand, composition, teloWindow=100, windowStep=6, plateauDetectionThreshold=-15, changeThreshold=-5, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=500, showGraphs=False):
    
    validate_seq_teloWindow(seq, teloWindow)

    if not isinstance(isGStrand, bool) and not isinstance(isGStrand, np.bool_):
        raise ValueError("isGStrand should be a boolean, or numpy boolean")

    if not isinstance(windowStep, int) or windowStep < 1:
        raise ValueError(
            "windowStep should be an int greater than or equal to 1")

    if not isinstance(targetPatternIndex, int) or targetPatternIndex > len(composition):
        raise ValueError(
            "targetPatternIndex should be an int and within the range of the composition list")

    if not isinstance(nucleotideGraphAreaWindowSize, int) or nucleotideGraphAreaWindowSize < 1:
        raise ValueError(
            "nucleotideGraphAreaWindowSize should be an int greater than or equal to 1")

    if not isinstance(showGraphs, bool):
        raise ValueError("showGraphs should be a boolean")

def validate_seq_teloWindow(seq, teloWindow):
    if not isinstance(teloWindow, int) or teloWindow < 6:
        raise ValueError(
            "teloWindow should be an int greater than or equal to 6")

    if len(seq) < teloWindow:
        raise Warning("sequence length is less than teloWindow")
