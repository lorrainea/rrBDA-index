rrBDA-index int
===

How to use
----------

### Installation

```
cd rrBDA-index_int
./pre-install.sh
make -f Makefile.32-bit.gc
```

### Usage

```
./rrbda-index_int <text_file> <ell> <pattern_file> <output_filename> <block_size> <index_filename>

<text_file> - name of input text file.
<ell> - lower bound on the length of input patterns to consider. 
<pattern_file> - name of input file containing the patterns.
<output_filename> - name of output file, where pattern occurrences will be output.
<block_size> - size of block to use for constructing the bd-anchors (bytes).
<index_filename> - name of the index file to be used (if it exists) otherwise to be created.
```

### Examples

```
 $ ./rrbda-index_ext ./data/text 3 ./data/patterns out 10 index
```
