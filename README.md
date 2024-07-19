rrBDA-index: Text Indexing for Long Patterns
===

Description
-----------
This repository maintains a time- and space-efficient construction algorithm of the <b>rrBDA-index</b>, a new variant of the <b>BDA-index</b> text index for long patterns introduced by [Loukides, Pissis, and Sweering](https://doi.org/10.1109/TKDE.2022.3231780).
This new construction relies on a linear-time algorithm for computing a <b>randomized</b> counterpart of <b>reduced</b> bd-anchors. 
It comes in two flavours: a semi-external memory implementation and an internal memory implementation.

Requirements
-----------
* A GNU/Linux system
* A modern C++17 ready compiler 

How to use
----------
<b>INPUT</b>: A file containing the text and a file containing a set of patterns seperated by a new line.

<b>OUTPUT</b>: A file containing the set of patterns and the starting position of their occurrences within the text.

### Installation

For the construction that uses internal memory only, compile as follows:
```
cd rrBDA-index_int
./pre-install.sh
make -f Makefile.32-bit.gcc
```

For the construction that uses both external and internal memory, compile as follows:
```
cd rrBDA-index_ext
./pre-install.sh
make -f Makefile.32-bit.gcc
```

### Usage

```
./rrbda-index_int <text_file> <ell> <pattern_file> <block_size> <output_filename> <index_filename>
./rrbda-index_ext <text_file> <ell> <pattern_file> <block_size> <ram_use> <output_filename> <index_filename>

<text_file> - name of input text file.
<ell> - lower bound on the length of input patterns to consider. 
<pattern_file> - name of input file containing the patterns.
<block_size> - size of block to use for constructing the bd-anchors (bytes).
<ram_use> - RAM usage for external SA and LCP array construction (MiB).
<output_filename> - name of output file, where pattern occurrences will be output.
<index_filename> - name of the index file to be used (if it exists) otherwise to be created.
```

### Examples

For the construction that uses internal memory only, run as follows:
```
 $ ./rrbda-index_int ./data/text 3 ./data/patterns 10 out index
```

For the construction that uses both external and internal memory, run as follows:
```
 $ ./rrbda-index_ext ./data/text 3 ./data/patterns 10 1024 out index
```
Note that the <b>exact same index</b> is constructed in the end.

Evaluation
-----------
We have conducted an extensive evaluation of different text indexes. We give a teaser below using the full human genome (version hg38) and 500K patterns (randomly selected from the text).

This plot depicts the size of the indexes for growing <b>ell</b>.
  <p align="center">
    <img src="https://github.com/lorrainea/rrBDA-index/blob/main/.images/hg38_size.png" alt="size" width=500 height=400>
  </p>

This plot depicts the average pattern matching time of the indexes for growing pattern length <b>|P|</b>.
  <p align="center">
    <img src="https://github.com/lorrainea/rrBDA-index/blob/main/.images/hg38_pattern.png" alt="query_time" width=500 height=400>
  </p>

Citation
--------

Lorraine A. K. Ayad, Grigorios Loukides, and Solon P. Pissis. 2024. Text indexing for long patterns using locally consistent anchors. [[arXiv]](https://arxiv.org/abs/2407.11819) (Under Review).

License
--------

GNU GPLv3 License; Copyright (C) 2024 Lorraine A. K. Ayad, Grigorios Loukides and Solon P. Pissis.
