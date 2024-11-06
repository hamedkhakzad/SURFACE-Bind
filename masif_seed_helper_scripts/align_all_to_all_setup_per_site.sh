#!/bin/bash

UND="_"

for d in out_peptides/*/*/; do
  echo $d
  dirname=$(basename $d)
  mkdir -p analysis_fixed_size$UND$dirname/
  mkdir -p analysis_fixed_size$UND$dirname/out_pdb/
  mkdir -p analysis_fixed_size$UND$dirname/out_data/
  find $d -name "*.score" -printf "%f\n" | sed 's/.score//g' > analysis_fixed_size$UND$dirname/list_out_peptides$UND$dirname.txt
  num_lines=$(wc -l analysis_fixed_size$UND$dirname/list_out_peptides$UND$dirname.txt | cut -d" " -f1)
  cp align_all_to_all_fixed_length.slurm align_all_to_all_fixed_length$UND$dirname.slurm
  sed -i "s/#SBATCH --array=.*/#SBATCH --array=1-$num_lines/" align_all_to_all_fixed_length$UND$dirname.slurm
done

