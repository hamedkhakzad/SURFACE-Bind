#!/bin/bash

## This script is to be run on Jed (CPU cluster)
## This script is to use ProteinMPNN to design a provided design(s)

python /work/lpdi/users/asadek/scaffold_opt_project/scripts/src/run_sol_mpnn_opt.py --rs_models_dir {Path to folder containing rosetta models} \
                                                                                    --pmpnn_repo {Path to ProteinMPNN repo} \
                                                                                    --path_to_conda {Path to source conda.sh, e.g. "/home/asadek/miniconda3/etc/profile.d/conda.sh"} \
                                                                                    --path_to_env {Path to activate conda environment, e.g. "/home/asadek/miniconda3/envs/ProteinMPNN"} \
                                                                                    --binder_chain 'B' \
                                                                                    --seqs_to_generate 25 \
                                                                                    --optimisation_strategy design_interface


