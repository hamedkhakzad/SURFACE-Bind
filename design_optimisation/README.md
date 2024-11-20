## Execution scripts for running the oprimisation pipeline 
*The scriptts are setup to run the pipeline on CPU / GPU clusters using SLURM as task manager*

**Software to install:** 
- ProteinMPNN (https://github.com/dauparas/ProteinMPNN.git)
- ColabFold   (https://github.com/sokrypton/ColabFold.git)
- ColabDesign (https://github.com/sokrypton/ColabDesign.git)

**Packages to install:**
- Biopython
- numpy
- pandas
- scipy
- matplotlib
- seaborn

------------------------------------------------------------------------------------------------------------------------------------------
## Run ProteinMPNN optimisation 
*provided example scripts for each setting: 01-scaffold optimisation while fixing the interacting residues, and 02-full binder design*

**Flags:**
rs_models_dir         : Path to folder containing rosetta designed models that will be input for optimisation. (str)
pmpnn_repo            : Path to ProteinMPNN repository. (str)
scripts_src_repo      : Path to scripts "src" folder.   (str)
binder_chain          : Binder chain ID to optimise.    (str)
interface_threshold   : Distance threshold for a residue to be considered a hotspot, Default 3.5Å. (float)
seqs_to_generate      : Number of sequences to generate per design. (int)
path_to_conda         : Path to source "conda.sh" file.     (str)
path_to_env           : Path to activate conda environment. (str)
optimisation_strategy : Name of ProteinMPNN strategy to be used, choices: ['fixed_interface', 'design_interface']
