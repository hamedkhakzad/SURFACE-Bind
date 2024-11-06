#!/bin/bash

addr="/work/lpdi/users/khakzad/masif_seed/masif_seed_search/data/masif_targets/targets/FINISHED_jobs"

cat list_jobs.txt | while read i; 
  do
    cd $addr/$i;
    find . -type f -name '*.score' -exec cat {} + > output_$i.score
    sort -nr -t: -k4 output_$i.score > output_sorted_$i.score
  done
