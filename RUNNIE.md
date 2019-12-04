![Oxford Nanopore Technologies logo](images/ONT_logo_590x106.png)


# Run-Length Encoded Basecaller (Runnie)
Runnie is an experimental basecaller that works in 'run-length encoded' space.
Rather than calling a sequence of bases, where one call is one base, runs of
bases (homopolymers) are called and so one call may represent many bases of the
same type.

Run-length encoding separates a sequence into two parts: a run-length
compressed sequencing containing only the identities of each run of bases, and
the corresponding length of each run.  Factorising in this manner has a
particular advantage for Oxford Nanopore reads since the identity of a
homopolymer is separated from its length.  Since homopolymer indels are a
significant class of error in single-read Nanopore basecalls, calls in "RLC" space
are more accurate and so may led to better mapping and consensus results.

    Sequence                CACCTTTTTTTGTAACGCTAAAGTCTCTTTTCAAACTTGCATTTTTGTAA
    Run-Length Compressed   CAC T      GTA CGCTA  GTCTCT   CA  CT GCAT    GTA
    Run-Lengths             112 7      112 11113  111114   13  12 1115    112


##  Installation
Assuming Flappie already works, see [README.md](README.md), then a `runnie` executable
can be created by:
```bash
make runnie
```

##  Usage
Basecalling using `runnie` is a two-stage process: firstly, a `.run` file is
created containing the Run-Length Compressed (RLC) basecall, along information
about the likely length of each run, then this is decoded into a **fasta**
basecall.  The executable `runnie` performs the initial RLC basecall, which is
then processed by the `misc/decode_runnie.py` script.

The `misc/decode_runnie.py` script requires a [numpy](https://www.numpy.org) installation.

```bash
#  ! It is highly recommended that OpenBLAS is run in single threaded mode
export OPENBLAS_NUM_THREADS=1
#  Process reads directory
runnie reads/ > basecalls.run
#  Create a .fasta from .run
python -O misc/decode_runnie.py basecalls.run > basecalls.fa 
#  Run in parallel
find reads -name \*.fast5 | parallel -P $(nproc) -X runnie | gzip > basecalls.run.gz
#  All in one go
find reads -name \*.fast5 | parallel -P $(nproc) -X runnie | python -O misc/decode_runnie.py --threads 4 > basecalls.fa
```

##  Run distributions
Rather than making a single call for the length of a homopolymer, `Runnie`
emits two parameters for each base called ('shape' and 'scale') that describe a
distribution over all possible lengths.  The estimate of the homopolymer
run-length is created from the distribution, e.g. by the mode or mean.

The advantage of describing the possible run-lengths by a distribution is that
it contains information about confidence in a particular run-length and likely
alternative calls.  The extra informartion provided by the distribution may be
particularly helpful for consensus calling, allowing different run-lengths to
be combined in a principled manner.

###  Discrete Weibull distribution
Currently, the distribution over run-lengths is assumed to be a [Discrete
Weibull
distribution](https://en.wikipedia.org/wiki/Discrete_Weibull_distribution) ,
with _shape_ and _scale_ parameters produced for each base by the basecalling
network.  The discrete Weibull distribution is formed by integating the
(continuous) [Weibull
distribution](https://en.wikipedia.org/wiki/Discrete_Weibull_distribution)
between neighbouring integers, so the probability that run-length is *l* is
*W(l + 1 | a, b) - W(l | a, b)*, where *W(x | a, b )* is the cumulative density
function for the Weibull distribution with shape *a* and scale *b*.

###  Estimating run-length
The exact proceedure for estimating the length of a homopolymer run from the
run-length distribution is being refined.  The current approach is to estimate
the mode of the Discrete Weibull distribution by rounding down the mode of the
equivalently parameterised continuous Weibull distribution, for which there is
a closed form expression.


##  Run format
The output of `runnie` is a simple encoding of the run-length compressed
basecall and information about the length of the run distribution.  While
having many deficiencies, the `.run` format is simple and contains richer
information about the run distribution than could be put in the **fasta** or
**fastq** calls.  

The primary purpose of the `.run` format is to allow experimentation with
different methods of decoding the run length and may be replaced by a different
format in the future.  It is unlikely that the ONT production basecallers will
produce this format if / when Run-Length Encoding is integrated into our
commercial products.

  
Each entry starts with a comment containing the read ID, followed by one line
for each base in the run-length compressed basecall.  These lines contain the
following values, tab-separated: base called, the _shape_ and _scale_ of the
run-length distribution, and the number of blocks in that run.

##  Limitations
Being a initial research release, `Runnie` several limitations and is only
fit for research use. In particular:

* Only R9.4.1 reads from MinION are supported.
* Modified basecalling is not supported.
* Multi-read fast5 files are not supported.

These limitations may be addressed in future releases.


## Licence and Copyright
(c) 2018 Oxford Nanopore Technologies Ltd.

Runnie is distributed under the terms of the Oxford Nanopore Technologies, Ltd.
Public License, v. 1.0.  If a copy of the License was not distributed with this
file, You can obtain one at http://nanoporetech.com



The vectorised math functions used by Flappie
[src/sse_mathfun.h](src/sse_mathfun.h) are from
http://gruntthepeon.free.fr/ssemath/ and the original version of this file is
under the 'zlib' licence.  See the top of
[src/sse_mathfun.h](src/sse_mathfun.h) for details.


Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools.  Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests.  However much as we would like to
rectify every issue and piece of feedback users may have, the developers may
have limited resource for support of this software.  Research releases may be
unstable and subject to rapid iteration by Oxford Nanopore Technologies.
