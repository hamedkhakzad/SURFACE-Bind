#!/bin/bash

## This script is to be run on Jed (CPU cluster)
## This script is to use ProteinMPNN to design a provided design(s) excpet the interaction interface within the provided cutoff

python /work/lpdi/users/asadek/scaffold_opt_project/scripts/src/run_sol_mpnn_opt.py --rs_models_dir {Path to folder containing rosetta models} \
                                                                                    --pmpnn_repo {Path to ProteinMPNN repo} \
                                                                                    --path_to_conda {Path to source conda.sh} \
                                                                                    --path_to_env {Path to activate conda environment} \
                                                                                    --scripts_src_repo {Path to scripts src repo} \
                                                                                    --binder_chain 'B' \
                                                                                    --interface_threshold 3.5 \
                                                                                    --seqs_to_generate 25 \
                                                                                    --optimisation_strategy fixed_interface


