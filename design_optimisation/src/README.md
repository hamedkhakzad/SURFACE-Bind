# Source scripts for Running design optimisation using ProteinMPNN and AlphaFold2 on a computing cluster

  - find_hotspots.py:
    This script parses a PDB file and determine the binder interface residues within a set distance cutoff from the               target (Default is set to 3.5Å). This script is necessary for running scaffold optimisation using ProteinMPNN.
  
  - run_sol_mpnn_opt.py:
    This script runs ProteinMPNN on input structures in two setting: 1. optimise non-interacting resdiues (i.e.,                  scaffold), and 2. full binder design.

  - run_af_monomer.py:
    This script runs AlphaFold2 - within ColabFold - to predict the optimised binder sequences for in silico folding              assessment.
    
  - run_af_monomer_analysis.py:
    This script analyses the AlphaFold2 predictions and filter them based on pLDDT and CA RMSD to model.


