import gzip
import io
import numpy as np
import os
# import pandas as pd
import torch

from Bio import PDB
from Bio import *
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# from dmasif_nlp.data import PairData
from torch.utils.data import Dataset, Subset
from tqdm import tqdm
from typing import Iterable

import pyrosetta
from pyrosetta import rosetta

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN

ELE2NUM = {"C": 0, "H": 1, "O": 2, "N": 3, "S": 4, "SE": 5}

def get_residue_coordinates(residue):
    coords = {}
    for atom in residue.get_atoms():
        coords[atom.name] = atom.get_coord()
    if 'CA' in coords:
        return coords['CA']
    if 'C' in coords:
        return coords['C']
    if 'N' in coords:
        return coords['N']
    return np.array(list(coords.values())).mean(axis=0)


def load_structure_gz(path):
    with gzip.open(path, 'rt') as f:
        handler = io.StringIO(f.read())
    return load_structure(handler)


def load_structure(path):
    parser = PDB.PDBParser()
    structure = parser.get_structure("structure", path)
    residues = structure.get_residues()
    atoms = structure.get_atoms()

    residuenumber = []
    residue2coords = {}
    for pos, residue in enumerate(residues):
        residue_number = residue.get_id()[1]
        residue_coords = get_residue_coordinates(residue)
        residue2coords[residue_number] = residue_coords
        residuenumber.append(residue_number)

    atom_coords = []
    atom_types = []
    for atom in atoms:
        atom_coords.append(atom.get_coord())
        atom_types.append(ELE2NUM[atom.element])

    atom_coords = np.stack(atom_coords)
    atom_types_array = np.zeros((len(atom_types), len(ELE2NUM)))
    for i, t in enumerate(atom_types):
        atom_types_array[i, t] = 1.0

    return atom_coords, atom_types_array, residue2coords, residuenumber


class Protein():
    def __init__(self, acc, coords, embeddings, labels, pdb_path, bs_threshold, device):
        
        ## the threshold we consider for MaSIF-site predictions to be considered site or not.
        self.bs_threshold = bs_threshold
        self.acc = acc
        self.coords = coords
        self.embeddings = embeddings
        self.labels = labels
        self.pdb_path = pdb_path
    
    def get_sequence(self):
        """
        Class function to get the sequence from the input pdb and store them in a list.
        """
        sequence = []
        target_pdb = self.pdb_path

        p = PDBParser()
        structure = p.get_structure("PDB0", target_pdb)
        model = structure[0]

        for chain in model:
            for residue in chain:
                sequence.append(residue.resname)

        return sequence

    @staticmethod
    def calc_sasa_per_res(pose):
        rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
        rsd_hydrophobic_sasa = pyrosetta.rosetta.utility.vector1_double()
        rosetta.core.scoring.calc_per_res_hydrophobic_sasa(pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4) # the probe radius
        return rsd_sasa

    def get_sasa(self):
        """
        Class function to calculate the solvent accessibility surface area.
        """
        pyrosetta.init('-ignore_unrecognized_res true')
        sasa_results = []
        target_pdb = self.pdb_path
        
        pose = rosetta.core.pose.Pose()
        rosetta.core.import_pose.pose_from_file(pose, target_pdb)
        sasa_results = self.calc_sasa_per_res(pose)
        return sasa_results

    def point2residue(self):
        """
        This function return a vector with the same length of surface_coords or the number of points and contain
        the index of the closest residue to each point. These indices are based on the residue_coord vector.
        """
        target_pdb = self.pdb_path
        atom_coords, atom_types, residue2coords, residue_number = load_structure(target_pdb)

        residuecoords = torch.tensor(
                        list(residue2coords.values()),  # considering all residues
                        dtype=torch.float32,
                        )

        coord = self.coords
        surface_coords = coord[:, None, :]
        residue_coords = residuecoords[None, :, :]
        distances_allvsall = torch.norm(surface_coords - residue_coords, dim=-1)  # (N, M)
        distances_min = torch.argmin(distances_allvsall, dim=1)
        return distances_min, residuecoords, residue_number

    def clustering_labels(self):
        """
        Function for clustering the point clouds that are in the binding sites.
        Basically it return the number of binding sites.
        """
        # target_pdb = os.path.join(self.pdb_path, self.acc[idx] + ".pdb")
        coord = self.coords
        label = self.labels

        sparse_label = label
        sparse_label[sparse_label >= self.bs_threshold] = 1
        sparse_label[sparse_label < self.bs_threshold] = 0

        masked_coord = sparse_label * coord
        a = masked_coord
        b = a[a.sum(dim=1).abs() != 0]

        data = b.numpy()
        model = DBSCAN(eps=2.5, min_samples=2)
        model.fit_predict(data)
        pred = model.fit_predict(data)

        ## get rid of cluster with tag 0 to avoide interfering with rest of calculation
        final_labels = model.labels_ + 1

        no_binding_sites = len(set(final_labels))
        #print("number of cluster found: {}".format(len(set(final_labels))))
        #print('cluster for each point: ', final_labels)

        norm_label = (final_labels - min(final_labels)) / (max(final_labels) - min(final_labels))

        mask = masked_coord.abs().sum(dim=1).bool()
        Y = torch.where(mask == True)[0].long()
        X = torch.tensor(final_labels).long()
        mask_all = torch.zeros(len(mask)).long()
        mask_all[Y] = X

        return no_binding_sites, b, norm_label, mask_all
