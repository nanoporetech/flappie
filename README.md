![Oxford Nanopore Technologies logo](images/ONT_logo_590x106.png)


# Flappie

## Overview

Basecall Fast5 reads using _flip-flop_ basecalling.  

## Features

* Flip-flop basecalling for the MinION platform
  * R9.4.1 (Native or PCR libraries)
  * R10C (PCR libraries only)
* Basecalling of 5mC in CpG context for R9.4.1

# Getting Started

## Input and Output

## Installation
Flappie has been tested on Ubuntu 16.04.5 LTS.  Other systems may be compatible.

Flappie models and other large resources are stored using [git lfs](https://git-lfs.github.com/) and this extension must be install to successfully clone the repository.

```bash
git clone https://github.com/nanoporetech/flappie
cd flappie
make flappie
```

### Compilation From Source
Flappie has the following dependences
* [CUnit](http://cunit.sourceforge.net/) library for unit testing
* [HDF5](https://www.hdfgroup.org/) library
* [OpenBLAS](https://www.openblas.net/) library for linear algebra


On Debian based systems, the following packages are sufficient (tested Ubuntu 14.04 and 16.04)
* libcunit1
* libcunit1-dev
* libhdf5
* libhdf5-dev
* libopenblas-base
* libopenblas-dev


## Usage

```bash
#  ! It is highly recommended that OpenBLAS is run in single threaded mode
export OPENBLAS_NUM_THREADS=1
#  List available models
flappie --model help
#  Basecall reads directory
flappie reads/ > basecalls.fq
#  Basecall using a different model
flappie --model r941_5mC reads/ > basecalls.fq
#  Output to SAM
flappie --format sam reads/ > basecalls.sam
#  Output to BAM
flappie --format sam reads | samtools view -Sb - > basecalls.bam
#  Dump trace data
flappie --trace trace.hdf5 reads > basecalls.fq
#  Basecall in parallel
find reads -name \*.fast5 | parallel -P $(nproc) -X flappie > basecalls.fq
#  Dump trace in parallel.  One trace per parallel process.
find reads -name \*.fast5 | parallel -P $(nproc) -X flappie --trace trace_{%}.hdf5 {} > basecalls.fq
```

# Help

## Licence and Copyright
(c) 2018 Oxford Nanopore Technologies Ltd.

Flappie is distributed under the terms of the Oxford Nanopore Technologies Developer licence.


The vectorised math functions used by Flappie [src/sse_mathfun.h](src/sse_mathfun.h) are from
http://gruntthepeon.free.fr/ssemath/ and the original version of this file is
under the 'zlib' licence.  See the top of [src/sse_mathfun.h](src/sse_mathfun.h) for details.


## FAQs

###  Platform support
The models contained contained in Flappie are trained using data from the MinION platform.  
Use on other platforms is not supported, although they may generalise to reads from the 
GridION platform due to the similarity of the hardware.
###  Methylation and other modifications
Flappie currently only calls 5mC methylation in CpG contexts.  Calling other modifications and 5mC is not currently supported.  Methylated calls are currently represented as a 'Z' base in the output
### Trace file
The trace information is output as a block x state matrix, where the states are the flip (uppercase) and flop (lowercase) bases in the order ACGTacgt or ACGTZacgtz for methylated calls.  The probabilities for each state are normalised into the range 0..255 and then represented by an unsigned 8bit integer.  Due to rounding, the sum of encoded probabilities for each blcok may not equal 255.
### Quality scores
The quality currently produced for FASTQ and SAM output are derived directly from the probabilistic model output by the _Flappie_ model and have not been calibrated.

## Abbreviations

## References and Supporting Information

