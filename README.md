# ONTbarcoder

1. [Overview](#overview)
2. [Getting the software](#getting-the-software)
3. [Installation](#installation)
4. [Types of analysis](#types-of-analysis)
5. [Test files](#test-files)


## Overview

ONTbarcoder (available from https://github.com/asrivathsan/ONTbarcoder) allows you to call DNA barcodes fom Nanopore sequence data. It allows the user to ONTbarcoder is designed to handle dual tagged amplicon pools. In the Conventional mode, it accepts a fastq file obtained after base-calling nanopore data and a demultiplexing file that specifies the tag combination for every sample. The [original paper](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01141-x) describes in detail how such tagged amplicons can be obtained and the pipeline implemented in ONTbarcoder. In the real-time barcoding mode, described [here](https://doi.org/10.1101/2023.06.26.546538), it accepts either fastq or raw files to call barcodes real-time, given user-defined demultiplexing information. 

ONTbarcoder has been optimized for protein coding gene like COI. While there are ways to obtain consensus for length variable non coding genes, this has not been extensively tested. 

## Getting the software
Please check the releases in this repository: https://github.com/asrivathsan/ONTbarcoder/releases. The ONTBarcoder_manual.pdf details more detailed guidelines for troubleshooting in case of issues. It is installed by unzipping the folder with the version of the program that supports the operating system on your computer. The folders are available from releases in Github. For MacOS, an .app bundle has been created. Kindly download the bundle relevant to the OS version. For further notes on MacOS compatibility and permissions: see the [MacOS section](#macintosh). Currently the software has been tested in OSX 12.5 and 13.5.

## Installation
#### Windows

Current version: Unzip [ONTbarcoder2.1.5_win.zip](https://github.com/asrivathsan/ONTbarcoder/releases/download/2.1.3/ONTbarcoder2.1.5_win.zip) and double click on the executable. In default mode the output files are stored in the same directory as the executable. Windows has a character limit of path length as 256 characters. Therefore, please avoid placing the software deep inside the directory structure. Avoid spaces in folder names

Older version: Unzip [ONTBarcoder_0.1.9_exe.win-amd64-2.7.zip](https://github.com/asrivathsan/ONTbarcoder/releases/download/0.1.9/ONTBarcoder_0.1.9_exe.win-amd64-2.7.zip) and double click on the executable. In default mode the output files are stored in the same directory as the executable

#### Macintosh

Use the DMG version [ONTbarcoder2.1.3.dmg](https://github.com/asrivathsan/ONTbarcoder/releases/download/2.1.3/ONTbarcoder2.1.3.dmg). In some tests, it was noted that permissions can be an issue in Mac. Transferring the app to Applications folder and running worked smoothly. If however one faces permission issues please modify System Preferences > Security and Privacy > Full Disk Access.

For older version of ONTbarcoder: Use the DMGs for the version corresponding to your Mac. We recommend using the upgraded Mac (macOS 11) and corresponding dmg file [ONTbarcoder_0.1.9_OS11.dmg](https://github.com/asrivathsan/ONTbarcoder/releases/download/0.1.9/ONTbarcoder_0.1.9_OS11.dmg). We have also compiled for older versions; [ONTbarcoder_0.1.9_OSX10_13.dmg](https://github.com/asrivathsan/ONTbarcoder/releases/download/0.1.9/ONTbarcoder_0.1.9_OSX10_13.dmg) has been tested in macOS 10.13 and 10.12. We would suggesting using this in first pass for 10.14/15

In case DMG files give problems due to permission, you may use the standalone.tgz: [macOS 11](https://github.com/asrivathsan/ONTbarcoder/releases/download/0.1.9/ONTbarcoder_0.9.1_OSX11_standalone.tgz), [macOS 10.13](https://github.com/asrivathsan/ONTbarcoder/releases/download/0.1.9/ONTbarcoder_0.1.9_OS10.13_standalone.tgz) files. These will contain numerous executables and the software should be run using ONTbarcoder_multiprocessing.

#### Linux
Current version: Uncompress the tgz archive [ONTbarcoder_2.1.3_linux.tgz](https://github.com/asrivathsan/ONTbarcoder/releases/download/2.1.3/ONTbarcoder_2.1.3_linux.tgz)
Start the software from the terminal as by 
```
cd directory_containing_ONTbarcoderfiles
./ONTbarcoder_multiprocessing
```

Older version: Uncompress the tgz archive [ONTbarcoder_multiprocessing_0.1.9_linux.tgz](https://github.com/asrivathsan/ONTbarcoder/releases/download/0.1.9/ONTbarcoder_multiprocessing_0.1.9_linux.tgz).  
## Types of analysis
There are two types of analysis implemented in the software
1. Conventional barcoding: This mode allows you to analyse a completed sequencing run. It has iterative mode that allows you to recover barcodes that do not meet initial quality controls. Conventional barcoding is described in [this](https://doi.org/10.1186/s12915-021-01141-x) study  and the detailed manual is provided in the github repository as [Conventional_barcoding_manual.pdf](https://github.com/asrivathsan/ONTbarcoder/blob/main/Conventional_barcoding_manual.pdf)
2. Real-time barcoding: This mode allows you to monitor the sequencing as the sequencing goes on in either mk1b or mk1c. It gives you rapid overview of your sequencing run and barcodes can be obtained within minutes of sequencing. Realtime barcoding is described [here](https://doi.org/10.1101/2023.06.26.546538) and the detailed manual is provided in the github repository as [Realtime_barcoding_manual.pdf](https://github.com/asrivathsan/ONTbarcoder/blob/main/Realtime_barcoding_manual.pdf)
3. Compare barcode sets: This module allows you to compare different barcode sets. The detailed manual is provided in this repository as [Compare_barcodes_manual.pdf](https://github.com/asrivathsan/ONTbarcoder/blob/main/Compare_barcodes_manual.pdf)

### Test files

Please use the Dataset_accessibility_and_descriptions file in the repository for access to test datasets. We recommend use of datasets A,B,C for initial testing before running the larger dataset.
