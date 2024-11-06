import sys
sys.path.append('dmasif_nlp')

import os
import pandas as pd
import numpy as np
import torch

from Bio import PDB
from Bio import *
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

#from data import load_structure
from dmasif_nlp.data import PairData
from dmasif_nlp.data_iteration import extract_single, process
from dmasif_nlp.geometry_processing import atoms_to_points_normals
from dmasif_nlp.model import dMaSIF
from dmasif_nlp.Arguments import parser
from torch_geometric.data import DataLoader
from tqdm import tqdm


parser.add_argument("--data", type=str, required=True)
parser.add_argument("--structures", type=str, required=True)
parser.add_argument("--dmasif_model", type=str, required=True)

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

    residue2coords = {}
    for pos, residue in enumerate(residues):
        residue_number = residue.get_id()[1]
        residue_coords = get_residue_coordinates(residue)
        residue2coords[residue_number] = residue_coords

    atom_coords = []
    atom_types = []
    for atom in atoms:
        atom_coords.append(atom.get_coord())
        atom_types.append(ELE2NUM[atom.element])

    atom_coords = np.stack(atom_coords)
    atom_types_array = np.zeros((len(atom_types), len(ELE2NUM)))
    for i, t in enumerate(atom_types):
        atom_types_array[i, t] = 1.0

    return atom_coords, atom_types_array, residue2coords


def generate_embeddings(protein_pair, net, args):
    protein_batch, _ = process(args, protein_pair, net)
    protein = extract_single(protein_batch, 0)
    outputs = net(protein)

    coords = outputs['P1']['xyz'].detach()
    normals = outputs['P1']['normals'].detach()
    emb_1 = outputs['P1']['embedding_1'].detach()

    labels = torch.sigmoid(outputs['P1']['iface_preds']).detach()

    return coords, normals, emb_1, labels


def main(args):
    # dataset.csv contains the accession numbers to access pdb files.
    data = pd.read_csv(os.path.join(args.data, 'dataset_all.csv'), sep=";")
    
    # Initializing model
    print('Loading dMaSIF model')
    net = dMaSIF(args)
    net.load_state_dict(torch.load(args.dmasif_model, map_location=args.device)["model_state_dict"])
    net = net.to(args.device)

    # Loading PDB files
    print('Loading PDB files')
    acc2protein_data = {}
    acc2residue_coords = {}
    for acc in tqdm(data['acc'].values):
        try:
            structure_path = os.path.join(args.structures, f'{acc}.pdb')
            if not os.path.exists(structure_path):
                print(f'Path {structure_path} does not exist')
                continue

            atom_coords, atom_types, residue2coords = load_structure(structure_path)
            atom_coords = torch.tensor(atom_coords, dtype=torch.float32, device=args.device)
            atom_types = torch.tensor(atom_types, dtype=torch.float32, device=args.device)
            batch_atoms = torch.zeros(len(atom_coords), dtype=torch.int8, device=args.device)
            point_coords, point_normals, batch_points = atoms_to_points_normals(
                atom_coords,
                batch_atoms,
                atomtypes=atom_types,
                resolution=args.resolution,
                sup_sampling=args.sup_sampling,
                distance=args.distance,
            )
            protein_data = PairData(
                pdb_code_p1=acc,
                atom_coords_p1=torch.tensor(atom_coords, dtype=torch.float32),
                atom_types_p1=torch.tensor(atom_types, dtype=torch.float32),
                xyz_p1=point_coords,
                normals_p1=point_normals,
            ).to(args.device)

            acc2protein_data[acc] = protein_data
            acc2residue_coords[acc] = residue2coords

        except (KeyError, IndexError, ValueError):
            continue

    # Running MaSIF-site and storing the data
    print('Running dMaSIF-site on input structures')
    batch_vars = ["xyz_p1", "atom_coords_p1"]
    loader = DataLoader(list(acc2protein_data.values()), batch_size=1, follow_batch=batch_vars)
    
    train_acc = pd.read_csv(os.path.join(args.data, 'dataset_all.csv'), sep=";").acc.unique()
    train_dataset = []

    for protein_pair in tqdm(loader):
        protein_pair = protein_pair.to(args.device)
        coords, normals, emb, labels = generate_embeddings(protein_pair, net, args)

        acc = protein_pair.pdb_code_p1[0]
        try:
            residue2coords = acc2residue_coords[acc]
            residue_coords = torch.tensor(
                list(residue2coords.values()),  # considering all residues
                dtype=torch.float32,
                device=args.device
            )
        except Exception as e:
            print(f'Problem with {acc}: {e}')
            continue

        result = {
            'acc': acc,
            'coords': coords,
            'normals': normals,
            'embeddings': emb,
	    'label': labels
        }
 
        if acc in train_acc: 
            train_dataset.append(result)
        else:
            print(f'Unknown acc={acc}')
       
    train_path = os.path.join(args.data, 'dataset_all_surfaces.pt')  
    torch.save(train_dataset, train_path)
  
if __name__ == '__main__':
    main(args=parser.parse_args())

