# This script is to define the hotspot residues based on distance from the interface 
# The default distance is 3.5A and in this version the Hydrogen atoms were removed from the calculations  

#libraries
import os
import sys
import Bio
import glob
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial.distance import cdist

#variables
pdb_folder = sys.argv[1]
ch_seed    = sys.argv[2]
thresh     = float(sys.argv[3])

#function
def get_hotspots(pdb_folder:str, binder_chain:str, thresh:float=3.5) -> str:
    """
    This function is to return the string of hotpots for a certain structure
    :params:
        - pdb_folder   : the folder containing the pdb to be processed, ONLY 1 PDB/FOLDER
        - binder_chain : the binder chain ID for which the hotspot residues will be defined
        - thresh       : Distance threshold below which a residue would be considered as a hotspot, DEFAULT = 3.5 A 
    """
    #define the path to the complex pdb
    pdb =[f for f in glob.iglob(os.path.join(pdb_folder, '*.pdb'))][0]
    
    #load the pdb file
    struct = PDBParser(QUIET=True).get_structure(os.path.basename(pdb), pdb)
    
    #define the chains
    binder_chains = list(binder_chain)
    target_chains = [x.id for x in struct.get_chains() if x.id not in binder_chains]
    
    #binder's chain length
    binder_chain_length = 0
    
    #binder's first residue number
    binder_first_res_num = 0
    
    #hotspots string
    hotspots =""
    
    #get binder length and first residue number
    for chain in struct[0]:
        if chain.id == binder_chains[0]:
            binder_chain_length = binder_chain_length + len(chain)
            binder_first_res_num = binder_first_res_num + (int(list(chain.get_residues())[0].id[1]))
    
    #get the atom coords 
    target_atoms = np.array([atom.get_coord() for atom in struct.get_atoms() if atom.get_full_id()[2] in target_chains if 'H' not in atom.get_full_id()[4][0]])
    binder_atoms = np.array([atom.get_coord() for atom in struct.get_atoms() if atom.get_full_id()[2] in binder_chains if 'H' not in atom.get_full_id()[4][0]])
    
    #map atoms to residues
    binder_residues = np.array([atom.get_parent().id[1] for atom in struct.get_atoms() if atom.get_full_id()[2] in binder_chains if 'H' not in atom.get_full_id()[4][0]])
    
    #generate the distance matrix
    dists = cdist(target_atoms, binder_atoms)
    
    #return the binder atoms that are closest to target chains
    closest_binder_atoms_clac = np.argmin(dists, axis=1)
     
    #get binder atoms that are within the targeted threshold
    closest_binder_atoms = closest_binder_atoms_clac[dists[np.arange(len(dists)), closest_binder_atoms_clac] < thresh]
    
    #get the hotspot residues that correspond to the closest atoms
    hotspot_residues = sorted(list(set(binder_residues[closest_binder_atoms])))
    
    #generate the hotspots list
    for item in hotspot_residues:
        hotspots = hotspots + str((int(item)-binder_first_res_num)+1) + " "
    hotspots = hotspots.rstrip(" ")
    
    return hotspots

#code
print(get_hotspots(pdb_folder=pdb_folder, binder_chain=ch_seed, thresh=thresh))
