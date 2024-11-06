#!/bin/bash

addr="/work/lpdi/users/khakzad/masif_seed/masif_seed_search/data/masif_targets/targets"

cpus=`squeue -u khakzad -h -t pending,running -r | wc -l`
cat list.txt | while read i; 
  do
  while [ $cpus -lt 2000 ]
    do
    name=${i%_A_beta*};
    cd $addr/$name"_A_beta"
    sbatch run_peptides.slurm $name"_A";
    sleep 5m
    done
  done