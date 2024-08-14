# TeloBP

TeloBP is an algorithm for identifying the telomere boundary in a genomic sequence. It scans through DNA sequence and finds the point at which the telomeric pattern breaks, making the telomere boundary point.

[![DOI](https://zenodo.org/badge/574249235.svg)](https://zenodo.org/doi/10.5281/zenodo.10826386)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Setup & Prerequisites

To use TeloBP, clone the repository and run this command to install TeloBP and its required packages:

```
pip install -e .
```

Alternatively, you can install the following packages manually:

```
Biopython
matplotlib
pandas
numpy
pandarallel
```

While TeloBP was developed in Python 3.10, the packages in the requirements.txt file should work with versions up to Python 3.12.

## Usage

### TeloBP Command Line Script

In the Scripts directory, the teloBPCmd.py script runs the teloBP analysis on a set of fastq files.
The script can be run from the command line using the following command:

```
python3 teloBPCmd.py <dataDir> <outputDir> [--fileMode] [--teloNP] [-v] [--targetQnamesCSV <csv file>] [--save_graphs]
```

The script takes the following arguments:
dataDir: The path to the data directory containing the fastq files to be analyzed.
OR the path to a single fastq file to be analyzed. Set the --fileMode flag to indicate this.
outputDir: The path to the output directory where the results will be saved.
--fileMode: Flag indicating that a single file is being analyzed.
--teloNP: Flag indicating that the teloNP analysis should be run instead of the teloBP analysis.
More details on the difference between TeloBP and TeloNP can be found below.
--targetQnamesCSV: Path to a csv file containing the qnames of the reads to be analyzed. The analysis will only be run on these reads.
--save_graphs: Flag indicating to generate graphs during analysis and save them to a pdf file. They will be saved in the main directory.
The program will run in single threaded mode when --save_graphs is set to True.
-v: Flag to enable verbose output.

The script will output a .csv file containing read qnames and telomere length values for all the reads which passed basic filtering.

### Genome Trimming

The "trimGenome.py" file has parameters preset for removing the telomeres on a genome. The script can be run with the following command:

```
python trimGenome.py <genome_file> <output_file> [output bed file]
```

If an output file path is giving as a third argument, a bed file containing telomere boundary coordinates for the ORIGINAL genome file will be created. The trimmed genome will simply cut the genome at these coordinates.

**NOTE**: This algorithm assumes that each chromosome is its own read, and not split into q and p arms.

### Genome Bed File Creator

This script will create a bed file containing the telomere boundary coordinates for a given genome. The script can be run with the following command:

```
python teloBPBedGenome.py <reference genome> <output bed file>
```

### TeloBP Function Arguments

Demo code for using TeloBP is provided in the demo.ipynb notebook file. The TeloBP function takes the following arguments:

```
getTeloBoundary(seq, isGStrand = None,  compositionGStrand=[], compositionCStrand = [], teloWindow=100, windowStep=6, plateauDetectionThreshold=-60, changeThreshold=-20, targetPatternIndex=-1, nucleotideGraphAreaWindowSize=500, showGraphs=False, returnLastDiscontinuity=False)
```

And returns the distance between the telomere boundary and the end of the sequence, or in other words, the length of the telomere.

While TeloBP's default parameters can be used to calculate a telomere boundary using just a sequence, depending on the application, finetuning the following parameters may be needed to improve accuracy. In most cases, the only arguments that need to be changed are the following:

**seq**: The sequence to be analyzed

**isGStrand**: Whether the sequence is the G-strand or not. Is None by default, and will be determined using the composition lists (count of G vs C repeats) if not specified.

**compositionGStrand and compositionCStrand**: The expected nucleotide composition of the telomere. This is used to calculate the expected telomere pattern. default for GStrand is ["GGG", 3/6], and ["CCC", 3/6] for CStrand. This means that the expected telomere pattern is "GGG" and should be found in 3/6 of the telomere.
Depending on the quality of the sequence, this may need to be changed. For example, nanopore reads with many misscalls may need a composition of ["GGG|AAA", 3/6, 3] to account for the misscalls.

When using regular expressions for the patterns, a third argument is required, specifying the length of the pattern being searched for. For example, ["GGG|AAA", 3/6, 3] would search for either "GGG" or "AAA", and expects to see them in 3/6 nucleotides in the telomere. The third argument specifies that the pattern being searched for is 3 nucleotides long. Another example would be ["TTAGGG|TTTGGG", 6/6, 6].

### TeloNP: TeloBP for Nanopore Reads

TeloNP uses the same TeloBP algorithm, but is pre-optimized for nanopore reads basecalled with guppy.

The current parameters are optimized to account for misscalls in Guppy basecalling.

#### TeloNP recommended usage

```
getTeloNPBoundary(seq)
```

Where seq is the sequence to be analyzed.

Optionally, you can specify the isGStrand parameter, which is a boolean value specifying whether the sequence is the G-strand or not. This may be preferable if the sequence strand is known through some previous analysis like an alignment, as it will save time by not having to calculate the composition of the sequence.

## TeloBP Algorithm Description

TeloBP works by scanning through the sequence and finding the point at which the telomeric pattern breaks, marking the telomere boundary point. It does this by scanning through the sequence in a window of size teloWindow, and calculates how similar the sequence is to the expected telomere composition. This similarity is calculated by counting the number of times the expected telomeric pattern appears in the window, and dividing it by the number of nucleotides in the window. The window is then moved along the sequence by the windowStep value, and the similarity is calculated again. This is repeated until the end of the sequence is reached. Graphing the offset scores produces a graph that looks like this:

![CHM13_chr1q_Nucleotides_Offset_Graph](https://github.com/CWGreider/GreiderLab/assets/78556850/82859bbd-158a-4534-ac4a-85fa7caeeec2)

Now we need to mark the point at which the offset score changes. To make sure that the change is not a random spike, we take the area under the curve of our offset scores, using the nucleotideGraphAreaWindowSize value as the window size. This produces a graph that looks like this:

![CHM13_chr1q_area_under_nucleotide_offset_graph](https://github.com/CWGreider/GreiderLab/assets/78556850/ffc42c8c-4c03-4e18-b5ae-5075ec4e8f17)

Here, we can expect to see the graph spike around the telomere boundary point. We scan ahead till we get an area value below the plateauDetectionThreshold, then start measuring the difference between area values till the slope plateaus. This is the point at which the telomeric pattern breaks, and is the telomere boundary point.

If "returnLastDiscontinuity" is set to True, we will look backwards till we reach a point above the changeThreshold value, then scan ahead till we get an area value below the plateauDetectionThreshold. This is useful for very noisy reads, where the first discontinuity may be a false positive.
