#!/bin/bash

## This script is to be run on a CPU cluster
## This script will select the top folded designs from  AF2 monomer predictions based on pLDDT>=80 and CA RMSD to model <= 1.5

python {PATH_to_src_DIR}/run_af2_monomer_analysis.py --mpnn_folder {Path to folder containing MPNN run outputs, e.g., "./sol_mpnn_fixed_interface/"} \
                                                     --rs_models_dir {Path to folder containing rosetta models} \
                                                     --rs_model_binder_chain 'B' \
                                                     --af_model_binder_chain 'A' \
                                                     --pLDDT_thresh 80.0 \
                                                     --RMSD_thresh 1.5 
