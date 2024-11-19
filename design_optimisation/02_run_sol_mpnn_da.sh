#!/bin/bash

## This script is to be run on a CPU cluster with SLURM managment system
## This script is to use ProteinMPNN to fully design provided design(s)
## In this example, we are generating 25 sequences for our binder chain "B"

python {PATH_TO_src_DIR}/src/run_sol_mpnn_opt.py --rs_models_dir {Path to folder containing rosetta models} \
                                                 --pmpnn_repo {Path to ProteinMPNN repo} \
                                                 --path_to_conda {Path to source conda.sh} \
                                                 --path_to_env {Path to activate conda environment} \
                                                 --binder_chain {Binder Chain ID} \
                                                 --seqs_to_generate 25 \
                                                 --optimisation_strategy design_interface


