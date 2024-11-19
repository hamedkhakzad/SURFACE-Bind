#!/bin/bash

## This script is to be run on Jed (CPU cluster)
## This script will select the top folded designs from  AF2 monomer based on pLDDT>=80 and CA RMSD <= 1.5 from binder
## Then these designs will be submitted for AF2 multimer predictions

python /work/lpdi/users/asadek/scaffold_opt_project/scripts/src/run_af2_monomer_analysis.py --af_pred_folder /work/lpdi/users/asadek/scaffold_opt_project/targets/pdl1/original_designs_af_monomer/af_pred/ \
                                                                                            --rs_models_dir /work/lpdi/users/asadek/scaffold_opt_project/targets/pdl1/5jds_selected_all/ \
                                                                                            --rs_model_binder_chain 'B' \
                                                                                            --af_model_binder_chain 'A'