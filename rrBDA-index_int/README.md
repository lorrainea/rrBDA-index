rrBDA-index using internal memory
===

How to use
----------

### Installation

```
cd rrBDA-index_int
./pre-install.sh
make -f Makefile.32-bit.gcc
```

### Usage

```
./rrbda-index_int <text_file> <ell> <pattern_file> <block_size> <output_filename> <index_filename>

<text_file> - name of input text file.
<ell> - lower bound on the length of input patterns to consider. 
<pattern_file> - name of input file containing the patterns.
<block_size> - size of block to use for constructing the bd-anchors (bytes).
<output_filename> - name of output file, where pattern occurrences will be output.
<index_filename> - name of the index file to be used (if it exists) otherwise to be created.
```

### Examples

```
 $ ./rrbda-index_int ./data/text 3 ./data/patterns 10 out index
```
