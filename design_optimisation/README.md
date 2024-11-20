## Execution scripts for running the oprimisation pipeline 
*The scriptts are setup to run the pipeline on CPU / GPU clusters using SLURM as task manager*

**Software to install:** 
- ProteinMPNN (https://github.com/dauparas/ProteinMPNN.git)
- ColabFold   (https://github.com/sokrypton/ColabFold.git)
- ColabDesign (https://github.com/sokrypton/ColabDesign.git)

**Packages to install:**
- scipy
- numpyall_sol_mpnn_designs
- pandas
- seaborn
- Biopython
- matplotlib

-----------------------------------------------------------------------------------------------------------------------------------------
## Run ProteinMPNN optimisation 
*Provided example scripts for each setting: 01-scaffold optimisation while fixing the interacting residues, and 02-full binder design* \
*Run on CPU cluster with SLURM managment system*

**Flags:**
- rs_models_dir         : Path to folder containing rosetta designed models that will be input for optimisation. (str)
- pmpnn_repo            : Path to ProteinMPNN repository. (str)
- scripts_src_repo      : Path to scripts "src" folder.   (str)
- binder_chain          : Binder chain ID to optimise.    (str)
- interface_threshold   : Distance threshold for a residue to be considered a hotspot, Default 3.5Å. (float)
- seqs_to_generate      : Number of sequences to generate per design. (int)
- path_to_conda         : Path to source "conda.sh" file.     (str)
- path_to_env           : Path to activate conda environment. (str)
- optimisation_strategy : Name of ProteinMPNN strategy to be used, choices: ['fixed_interface', 'design_interface']

**output folders:**
- _temp : folder containing temporary files for parallization, emptied when run finish.
- slurms: folder containing slurm files submitted to the CPU cluster.
- sol_mpnn_{strategy_name}_input  :  folder containing sorted files that will be passed as input for ProteinMPNN, compressed after usage.
- sol_mpnn_{strategy_name}_output : folder containing the ProteinMPNN output.

**output files:**
- all_sol_mpnn_input_des_paths_design_interface.txt : file containing the paths for all ProteinMPNN inputs.

**recommendations for running the script:**
- Prepare the run folder and copy the excution scripts to it.
- The execution files are simply run using "source 01_run_sol_mpnn_fh.sh".
- Avoid "_" and "x" in naming the rosetta models as they are part of string splitting and identification throughout the pipeline. 

-----------------------------------------------------------------------------------------------------------------------------------------
## Run AlphaFold2 assessment of optimised sequences in-silico folding  
*This step predicts the design alone in "single sequence" mode*
*Provided example script: 03_run_af_monomer.sh* \
*Run on GPU cluster with SLURM managment system*

**Flags:**
- mpnn_folder   : Path to ProteinMPNN run output folder. (str)
- num_of_seqs   : Number of unique optimised sequences to predict per design selected based on ProteinMPNN global score.  (int)
- path_to_conda : Path to source "conda.sh" file.     (str)
- path_to_env   : Path to activate conda environment. (str)

**output folders:**
- sol_mpnn_{strategy_name}_output_analysis : Folder containing the ProteinMPNN output analysis.
- af_monomer_slurms : Folder containing slurm files submitted to the GPU cluster.
- af2_monomer_fa_files : Fxsolder containing FASTA files of optimised sequences to predict.
- af2_monomer_pred : folder containing AlphaFold2 predictions.

**output files:**
- all_sol_mpnn_designs.csv      : File containing ProteinMPNN scores for all optimised designs.
- selected_sol_mpnn_designs.csv : File containing ProteinMPNN scores for the selected optimised designs to predict with AlphaFold2.
- fasta_paths.list : File containing paths for FATSA files submitted for AlphaFold2 prediction.
- {DATE}_pmpnn_scores.png (.svg) : Plotting ProteinMPNN scores for the optimised sequences, Selected sequences (Coloured), All sequences (Gray).

