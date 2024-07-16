rrBDA-index ext
===

How to use
----------

### Installation

```
cd rrBDA-index_ext
./pre-install.sh
make -f Makefile.32-bit.gcc
```

### Usage

```
./rrbda-index_ext <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>

<text_file> - name of input text file.
<ell> - lower bound on the length of input patterns to consider. 
<pattern_file> - name of input file containing the patterns.
<output_filename> - name of output file, where pattern occurrences will be output.
<ram_use> - RAM usage for external SA and LCP array construction (MiB).
<block_size> - size of block to use for constructing the bd-anchors (bytes).
<index_filename> - name of the index file to be used (if it exists) otherwise to be created.
```

### Examples

```
 $ ./rrbda-index_ext ./data/text 3 ./data/patterns out 1024 10 index
```
