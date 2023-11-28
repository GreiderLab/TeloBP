# TeloBP

TeloBP is an algorithm for identifying the telomere boundary in a genomic sequence. It scans through DNA sequence and finds the point at which the telomeric pattern breaks, making the telomere boundary point.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

TeloBP requires the following packages to be installed:

```
Biopython
matplotlib
pandas
numpy
```

### Setup

To use TeloBP, clone the repository and run this command to install the required packages:

```
pip install -r requirements.txt
```

## Usage

### Genome Trimming

The scripts folder contains the script "trimGenome.py" which has parameters preset for removing the telomeres on a genome. The script can be run with the following command:

```
python trimGenome.py <genome_file> <output_file>
```

**NOTE**: This algorithm assumes that each chromosome is its own read, and not split into q and p arms.

### TeloBP Function Arguments

Demo code for using TeloBP is provided in the demo.ipynb notebook file. The TeloBP function takes the following arguments:

```
getTeloBoundary(seq, isGStrand, composition=[], teloWindow=100, windowStep=6, plateauDetectionThreshold=-60, changeThreshold=-20, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=500, showGraphs=False, returnLastDiscontinuity=False)
```

And returns the distance between the telomere boundary and the end of the sequence, or in other words, the length of the telomere.

In most cases, the only arguments that need to be changed are the following:

**seq**: The sequence to be analyzed

**isGStrand**: Whether the sequence is the G-strand or not

**composition**: The expected nucleotide composition of the telomere. This is used to calculate the expected telomere pattern. default: ["GGG", 3/6]. This means that the expected telomere pattern is "GGG" and should be found in 3/6 of the telomere.
Depending on the quality of the sequence, this may need to be changed. For example, nanopore reads with many misscalls may need a composition of ["GGG|AAA", 3/6, 3] to account for the misscalls.

When using regular expressions for the patterns, a third argument is required, specifying the length of the pattern being searched for. For example, ["GGG|AAA", 3/6, 3] would search for either "GGG" or "TTT", and expects to see them in 3/6 nucleotides in the telomere. The third argument specifies that the pattern being searched for is 3 nucleotides long.

### TeloNP: TeloBP for Nanopore Reads

TeloNP uses the same TeloBP algorithm, but is pre-optimized for nanopore reads.

The current parameters are optimized to account for misscalls in Guppy basecalling.

#### TeloNP recommended usage

```
getTeloNPBoundary(seq, isGStrand)
```

Where seq is the sequence to be analyzed, and isGStrand is a boolean value specifying whether the sequence is the G-strand or not.

## TeloBP Algorithm Description

TeloBP works by scanning through the sequence and finding the point at which the telomeric pattern breaks, marking the telomere boundary point. It does this by scanning through the sequence in a window of size teloWindow, and calculates how similar the sequence is to the expected telomere composition. This similarity is calculated by counting the number of times the expected telomeric pattern appears in the window, and dividing it by the number of nucleotides in the window. The window is then moved along the sequence by the windowStep value, and the similarity is calculated again. This is repeated until the end of the sequence is reached. Graphing the offset scores produces a graph that looks like this:

![CHM13_chr1q_Nucleotides_Offset_Graph](https://github.com/CWGreider/GreiderLab/assets/78556850/82859bbd-158a-4534-ac4a-85fa7caeeec2)

Now we need to mark the point at which the offset score changes. To make sure that the change is not a random spike, we take the area under the curve of our offset scores, using the nucleotideGraphAreaWindowSize value as the window size. This produces a graph that looks like this:

![CHM13_chr1q_area_under_nucleotide_offset_graph](https://github.com/CWGreider/GreiderLab/assets/78556850/ffc42c8c-4c03-4e18-b5ae-5075ec4e8f17)

Here, we can expect to see the graph spike around the telomere boundary point. We scan ahead till we get an area value below the plateauDetectionThreshold, then start measuring the difference between area values till the slope plateaus. This is the point at which the telomeric pattern breaks, and is the telomere boundary point.

If "returnLastDiscontinuity" is set to True, we will look backwards till we reach a point above the changeThreshold value, then scan ahead till we get an area value below the plateauDetectionThreshold. This is useful for very noisy reads, where the first discontinuity may be a false positive.
