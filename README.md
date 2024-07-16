rrBDA-index: Text Indexing for Long Patterns using Randomized Reduced Anchors
===


How to use
----------
<b>INPUT</b>: A file containing a single text and a file containing a set of patterns seperated by a new line.

<b>OUTPUT</b>: A file containing the set of patterns and the starting position of their occurrences within the text.

### Installation

```
cd rrBDA-index_int
./pre-install.sh
make -f Makefile.32-bit.gcc
```

```
cd rrBDA-index_ext
./pre-install.sh
make -f Makefile.32-bit.gcc
```

### Usage

```
./rrbda-index_int <text_file> <ell> <pattern_file> <output_filename> <block_size> <index_filename>
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
 $ ./rrbda-index_int ./data/text 3 ./data/patterns out 10 index
 $ ./rrbda-index_ext ./data/text 3 ./data/patterns out 1024 10 index
```
