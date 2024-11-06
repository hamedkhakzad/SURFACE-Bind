#!/bin/bash

addr="/work/lpdi/users/khakzad/masif_seed/masif_seed_search/data/masif_targets/targets/FINISHED_jobs"

UND="_"

cat list_jobs.txt | while read i; 
  do
    cd $i
    cp ../../O15399_A/align_all_to_all_fixed_length.slurm .
    cp ../../O15399_A/align_all_to_all_setup_per_site.sh .
    ./align_all_to_all_setup_per_site.sh
    pwd
    sleep 1m
    for d in $(cd $addr/$i/out_peptides/$i/ && echo */); do
      site=$(echo basename "$d" | awk -F'[_/]' '{print $2}')
      echo $site
      echo "sbatch align_all_to_all_fixed_length_site_"$site".slurm" $i $site
      sbatch align_all_to_all_fixed_length_site_$site.slurm $i $site
      echo "alignment is succesfully finished!"
      sleep 5m
      sbatch seed_clustering.slurm $site
      echo "clustering is succesfully finished!"
      sleep 1m
    done
    cd ../
  done
