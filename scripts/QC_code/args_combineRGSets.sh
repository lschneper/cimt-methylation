#!/bin/bash

# Loop through CVH RGSet subsets and combine them with the old data. Start at 2 (first one was previously combined )!
for index in $(seq 2 29);
do
	Rscript args_combineRGSets.R ${index}
done 
