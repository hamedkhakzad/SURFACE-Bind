#!/bin/bash

## This script is to be run on Izar (GPU cluster)
## This script will select the top 5 MPNN designs based on the global score and submit the results to AF2 for monomer prediction using single sequence mode

python /work/lpdi/users/asadek/scaffold_opt_project/scripts/src/run_af_monomer.py --mpnn_folder {Path to folder containing MPNN run outputs} \
                                                                                  --num_of_seqs 5 \
                                                                                  --path_to_conda {Path to source conda.sh, e.g. "/home/asadek/miniconda3/etc/profile.d/conda.sh"} \
                                                                                  --path_to_env {Path to activate conda environment, e.g. "/home/asadek/miniconda3/envs/ColabFold"} \
