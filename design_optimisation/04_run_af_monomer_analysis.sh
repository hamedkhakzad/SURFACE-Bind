#!/bin/bash

## This script is to be run on a CPU cluster
## This script will select the top folded designs from  AF2 monomer predictions based on pLDDT>=80 and CA RMSD to model <= 1.5

python /work/lpdi/users/asadek/scaffold_opt_project/scripts/src/run_af2_monomer_analysis.py --af_pred_folder /work/lpdi/users/asadek/scaffold_opt_project/targets/pdl1/original_designs_af_monomer/af_pred/ \
                                                                                            --rs_models_dir {Path to folder containing rosetta models} \
                                                                                            --rs_model_binder_chain 'B' \
                                                                                            --af_model_binder_chain 'A'
