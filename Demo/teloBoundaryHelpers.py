import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    for i in range(0, len(row)-windowSize, 1):
        area = row[i:i+windowSize].sum()
        areaList.append((area/windowSize))
        # areaList.append(area)
    return areaList
    

def graphLine(rowIn, labelIn, windowStep, boundaryPoint=-1):
    # smoothing_window = 12
    # smoothed_data = np.zeros_like(rowIn)
    # smoothed_data = np.convolve(rowIn, np.ones(smoothing_window), mode='same') / smoothing_window
    
    # derivatives = np.gradient(smoothed_data, axis=0)
    # row=derivatives
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
    # Make a plot with 4 lines
    # x axis is the index of the offset
    # y axis is the offset value
    # 4 lines for T, A, G, C

    # Example multidimensional array
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
        ax.plot(x, transposed_data[i, :], label=f"{labelDict[i]} Offset")

    # Set labels and title
    ax.set_xlabel('Distance from end of sequence')
    ax.set_ylabel('Nucleotide offset')
    ax.set_title('Nucleotide Offsets from Expected Telomere Composition')
    ax.legend()
    plt.show()



manualLabels = {
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

def testTeloGenomePosition(chr, pos, testDict = manualLabels):
    # print(chr)
    if chr not in testDict:
        print("Error: chromosome not found in dictionary, returning 0")
        return 0
    expectedPos = testDict[chr]
    # print("Expected position: " + str(expectedPos))
    print(chr + " offset: " + str(pos - expectedPos)+"bp (obs - exp)")
    return (pos - expectedPos)

def testTeloLength(chr, length, testDict = manualLabels):
    if chr not in testDict:
        print("Error: chromosome not found in dictionary, returning 0")
        return 0
    expectedLength = testDict[chr]
    print(chr + " offset: " + str(length - expectedLength)+"bp (obs - exp)")
    return (length - expectedLength)


def write_bed_file(file_path, bed_data):
    with open(file_path, 'w') as bed_file:
        for entry in bed_data:
            bed_file.write('\t'.join(str(e) for e in entry) + '\n')

def validate_parameters(seq, isGStrand, composition, teloWindow = 100, windowStep = 6, maxAreaThreshold = -15, minAreaThreshold = -5, targetPatternIndex=-1, nucleotideGraphAreaWindowSize = 500, showGraphs = False):
    if not isinstance(teloWindow, int) or teloWindow < 6:
        raise ValueError("teloWindow should be an int greater than or equal to 6")
    
    if len(seq) < teloWindow:
        raise ValueError("Error: sequence length must be greater than teloWindow")
    
    if not isinstance(isGStrand, bool):
        raise ValueError("isGStrand should be a boolean")

    if not isinstance(windowStep, int) or windowStep < 1:
        raise ValueError("windowStep should be an int greater than or equal to 1")
    
    if not isinstance(targetPatternIndex, int) or targetPatternIndex > len(composition):
        raise ValueError("targetPatternIndex should be an int and within the range of the composition list")
    
    if not isinstance(nucleotideGraphAreaWindowSize, int) or nucleotideGraphAreaWindowSize < 1:
        raise ValueError("nucleotideGraphAreaWindowSize should be an int greater than or equal to 1")
    
    if not isinstance(showGraphs, bool):
        raise ValueError("showGraphs should be a boolean")
    

