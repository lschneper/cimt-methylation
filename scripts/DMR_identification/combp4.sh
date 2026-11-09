#!/bin/bash

# For Manhattan step, use filter_and_annotate_rgns.Rmd first!


comb-p acf -d 1:1500:200 -c 4 $1.bed > tmp.acf.1500-200.txt
comb-p slk --acf tmp.acf.1500-200.txt -c 4 $1.bed > tmp.slk.1500-200.txt
comb-p fdr tmp.slk.1500-200.txt > tmp.fdr.1500-200.txt
comb-p peaks --dist 1500 --seed 0.1 tmp.fdr.1500-200.txt > tmp.peak.1500-200.txt
comb-p region_p -p $1.bed -r tmp.peak.1500-200.txt -s 200 -c 4 > $1.regions.sig.1500-200.bed
# comb-p manhattan -c 4 --regions $1.regions.sig.1500-200.top_ann.bed $1.bed --image $1.manhattan.1500-200.png
