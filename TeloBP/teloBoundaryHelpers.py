import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

from constants import manualLabelsCHM13Positions, manualLabelsCHM13


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


def graphLine(rowIn, labelIn, windowStep, boundaryPoint=-1):
    row = rowIn
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
        raise ValueError(
            "Error: sequence length must be greater than teloWindow")
