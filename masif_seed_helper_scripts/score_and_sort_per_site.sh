#!/bin/bash

addr="/work/lpdi/users/khakzad/masif_seed/masif_seed_search/data/masif_targets/targets/FINISHED_jobs"
UND="_"

cat list_jobs.txt | while read i; 
  do
    name=${i%_A*};
    for d in $addr/$i/out_peptides/$name"_A"/*/; do
      cd $d
      result=${PWD##*/}
      find . -type f -name '*.score' -exec cat {} + > output_$i$UND$result.score
      sort -nr -t: -k4 output_$i$UND$result.score > output_sorted_$i$UND$result.score
    done
  done
