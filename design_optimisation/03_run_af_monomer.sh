#!/bin/bash

## This script is to be run on GPU cluster with with SLURM managment system
## This script will select the top 5 unique ProteinMPNN sequences/design based on the global score and submit the results to AF2 for monomer prediction using single sequence mode
## This script uses ColabFold for running the predictions

python {PATH_to_src_DIR}/src/run_af_monomer.py --mpnn_folder {Path to folder containing MPNN run outputs, e.g., "./sol_mpnn_fixed_interface/"} \
                                               --num_of_seqs 5 \
                                               --path_to_conda {Path to source conda.sh} \
                                               --path_to_env {Path to activate conda environment} \
