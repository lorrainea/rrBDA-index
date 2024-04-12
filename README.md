rBDA-index_I_int, rBDA-index_II_int, rBDA-index_I_ext and rBDA-index_II_ext
===

<b>Installation</b>: To install and compile BDA-index_I or BDA-index_II, read the INSTALL file within the BDA-index_I_int, BDA-index_II_int, BDA-index_I_ext or BDA-index_II_ext folders.

<b>INPUT</b>: A file containing a single text and a file containing a set of patterns seperated by a new line.

<b>OUTPUT</b>: A file containing the set of patterns and the starting position of their occurrences within the text.


```
Usage: 
rBDA-index_I_int
./rbda-index_I <text_file> <ell> <pattern_file> <block_size> <output_filename> <index_filename>

rBDA-index_II_int
./rbda-index_II <text_file> <ell> <pattern_file> <block_size> <output_filename> <index_filename>

rBDA-index_I_ext
./rbda-index_I <text_file> <ell> <pattern_file> <block_size> <ram_use> <output_filename> <index_filename>

rBDA-index_II_ext
./rbda-index_II <text_file> <ell> <pattern_file> <block_size> <ram_use> <output_filename> <index_filename>


<text_file> - name of input text file
<ell> - minimum size of pattern to consider searching for within text. 
<pattern_file> - name of input file containing patterns
<output_filename> - name of output file where pattern occurrences will be placed.
<ram_use> - ram usage for external SA and LCP
<block_size> - size of block size b to use.
<index_filename> - name of the index file to be used if they exist otherwise to be created.
```

<b>Examples</b>
```
rBDA-index_I_int
 $ ./rbda-index_I ./data/text 3 ./data/patterns 10 out 150 index

rBDA-index_II_int
 $ ./rbda-index_II ./data/text 3 ./data/patterns 10 out 150 index
 
rBDA-index_I_ext
 $ ./rbda-index_I ./data/text 3 ./data/patterns 10 2000 out 150 index

rBDA-index_II_ext
 $ ./rbda-index_II ./data/text 3 ./data/patterns 10 2000 out 150 index
```
