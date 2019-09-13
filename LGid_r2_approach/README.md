# LD-based approach to identify variant position on the linkage map

One script that does it all (I know, this is often a bad sign but we are working hard to solve this)!  
```
Rscript scripts/LD_pure_inds.R -V data/genotype.table -D data -E wave -I CZA -L 6
```

Unfortunately, the file `data/genotype.table` is missing but Samuel can provide that in no time. Otherwise, if you are too eager and cannot wait for it, you can feed the `-V` argument with any file that looks like this:

| CHROM | POS | TYPE | REF | ALT | CZA001.GT | CZA002.GT | CZA100.GT | ... |
| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| Contig100005 | 1204 | SNP | G | A | G/A | A/A | G/A | ... |
| Contig99939 | 3754 | INDEL | GTCA | G | G/G | GTCA/GTCA | GTCA/GTCA | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |

The current version of the script works only for the island CZA and this is why the directory `data` contains only files recorded for CZA. If somebody would like to run some tests on other islands, this is what he/she/it must do:  
1. Retrieve data for the other islands. The filenames should follow the same pattern as the one for CZA.  
1. Comment out lines from 147 to 149 of `scripts/LD_pure_inds.R`.  
1. Comment in lines from 150 to 169 and from 173 to 179 of `scripts/LD_pure_inds.R`.

GOOD LUCK!
