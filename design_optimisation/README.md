## Execution scripts for running the oprimisation pipeline 
*The scriptts are setup to run the pipeline on CPU / GPU clusters using SLURM as task manager*

**Software to install:** 
- ProteinMPNN (https://github.com/dauparas/ProteinMPNN.git)
- ColabFold   (https://github.com/sokrypton/ColabFold.git)
- ColabDesign (https://github.com/sokrypton/ColabDesign.git)

**Packages to install:**
- scipy
- numpy
- pandas
- seaborn
- Biopython
- matplotlib

-----------------------------------------------------------------------------------------------------------------------------------------
## I. Run ProteinMPNN optimisation 
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

**remarks:** \
These parameters were kept, for modifying them, you need to alter them in the SLURM file writing part inside "run_sol_mpnn_opt.py" script \
" --sampling_temp "0.1"  --seed 37  --batch_size 1 "

**recommendations for running the script:**
- Prepare the run folder and copy the excution scripts to it.
- The execution files are simply run using "source 01_run_sol_mpnn_fh.sh".
- Avoid "_" and "x" in naming the rosetta models as they are part of string splitting and identification throughout the pipeline. 

-----------------------------------------------------------------------------------------------------------------------------------------
## II. Run AlphaFold2 assessment of optimised sequences in-silico folding  
*This step predicts the design sequence only in "single sequence" mode* \
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

-----------------------------------------------------------------------------------------------------------------------------------------
## III. Filtering AlphaFold2 predictions
*This step filters the AlphaFold2 predictions based on pLDDT (Default >= 80) and C⍺ RMSD to model (Default <= 1.5)* \
*Provided example script: 04_run_af_monomer_analysis.sh* \
*Run on CPU cluster, inteactive*

**Flags:**
- mpnn_folder           : Path to ProteinMPNN run output folder. (str)
- rs_models_dir         : Path to folder containing rosetta designed models. (str)
- rs_model_binder_chain : Rosetta model binder chain ID. (str)
- af_model_binder_chain : AlphaFold2 model chain ID. (str)
- pLDDT_thresh          : pLDDT threshold for filtering, Default 80.
- RMSD_thresh           : C⍺ RMSD to model threshold for filtering, Default 1.5Å.

**output folders:**
- af2_monomer_pred_analysis : Folder containing the AlphaFold2 predictions analysis.
- aligned_models            : Folder containing aligned AlphaFold2 to Rosetta binder models.
- sel_designs               : Folder containing the filter passing AlphaFold2 models.

**output files:**
- af2_monomer_all_decoys_analyzed.csv          : File containing the scores for all AlphaFold2 predictions.
- af2_monomer_all_analyzed_decoys_averaged.csv : File containing the averaged scores for all AlphaFold2 predictions. Here we averaged the scores of the 5 AlphaFold2 models for a given sequence.
- af2_monomer_best_decoys.csv                  : File containig the scores of the AlphaFold2 best models, the ones passing the pLDDT threshold.
- af2_monomer_folded_decoys_rmsd.csv           : File containing the C⍺ RMSD to model values for the AlphaFold2 best models.
- af2_monomer_filtered_decoys_by_rmsd.csv      : File containing the C⍺ RMSD to model values for the AlphaFold2 models passing the C⍺ RMSD to model threshold.
- af_monomer_selection.csv                     : File containing the ProteinMPNN, AlphaFold2 and sequence properties of passing designs.

-----------------------------------------------------------------------------------------------------------------------------------------
## IV. Interface quality assessment using AlphaFold2
*This step filters the selected optimised designs based on ipTM (Default >= 0.7) \
*Provided example notebook: 05_reperdict_sequences_with_AF2_and_parital_templates.ipynb* \
*Run on GPU cluster, inteactive*

Please follow instructions in the notebook. 
