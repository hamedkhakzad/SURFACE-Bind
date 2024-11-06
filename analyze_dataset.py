import argparse
import os
from data import Protein
from tqdm import tqdm
import pandas as pd
import nglview as ng
from nglview.color import ColormakerRegistry
from pdbparser.pdbparser import pdbparser
from sklearn import decomposition
import torch

import sys
import numpy as np


AAKD2NUM = {"ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
            "GLU": -3.5, "GLN": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
            "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
            "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2}


np.set_printoptions(threshold=sys.maxsize)

parser = argparse.ArgumentParser()
parser.add_argument("--data", type=str, required=True)
parser.add_argument("--structures", type=str, required=True)
parser.add_argument("--outpath", type=str, required=True)
parser.add_argument("--threshold", type=float, required=False)
parser.add_argument("--noPointscluster", type=int, required=False)
parser.add_argument("--device", type=str, required=True)

def flat_score(coords_xyz):
    pca = decomposition.PCA(n_components=3)
    pca.fit(coords_xyz)
    pca_r = pca.explained_variance_ratio_
    
    return pca_r

def save_pointcloud(main_pdb, acc, coord_file, labels_file, out_path):
    # Normalize embedding to represent a b-factor value between 0-100
    b_factor = []
    for emb in labels_file:
        b_factor.append(emb)
    
    # writing a psudo pdb of all points using their coordinates and H atom.
    records = []

    for i in range(len(coord_file)):
        points = coord_file[i]
        x_coord = points[0]
        y_coord = points[1]
        z_coord = points[2]

        records.append( { "record_name"       : 'ATOM',
                    "serial_number"     : len(records)+1,
                    "atom_name"         : 'H',
                    "location_indicator": '',
                    "residue_name"      : 'XYZ',
                    "chain_identifier"  : '',
                    "sequence_number"   : len(records)+1,
                    "code_of_insertion" : '',
                    "coordinates_x"     : x_coord,
                    "coordinates_y"     : y_coord,
                    "coordinates_z"     : z_coord,
                    "occupancy"         : 1.0,
                    "temperature_factor": b_factor[i]*100,
                    "segment_identifier": '',
                    "element_symbol"    : 'H',
                    "charge"            : '',
                    } )

    pdb = pdbparser()
    pdb.records = records

    pdb.export_pdb(os.path.join(out_path,acc+"_pointcloud_clusters.pdb"))
    return 1

def main(args):

    # loading data from the dataset that prepared before and contain MaSIF-site predictions.
    dataset_path = os.path.join(args.data, 'dataset_surfaces.pt')
    pdb_path = os.path.join(args.structures)
    if os.path.exists(dataset_path):
        data = torch.load(dataset_path, map_location=torch.device(args.device))
    else:
        print("No such file or directory", dataset_path)

    # we are going to store results in this array and later write them as csv file. Each entry will be a dictionary.
    metrics = []
    # metrics.append({'pdb': pdb, 'roc_auc': roc_auc})

    for item in data:
	# initializing all the outputs that we are going to produce per protein and store in the dictionary
        sequence = []
        sasa = []
        point2residue = []
        residue2coords = []
        residue_number = []
        no_binding_sites = 0
        filtered_site_points = []
        normalized_labels = []
        cluster_for_all_points = []

        acc = item['acc']
        coords = item['coords']
        embeddings = item['embeddings']
        labels = item['label']

        target_pdb = os.path.join(pdb_path,acc+".pdb")
	
        protein = Protein(acc, coords, embeddings, labels, target_pdb, args.threshold, args.device)

        sequence = protein.get_sequence()
        sasa = protein.get_sasa()
        point2residue, residue2coords, residue_number = protein.point2residue()
        no_binding_sites, filtered_site_points, normalized_labels, cluster_for_all_points = protein.clustering_labels()

        ## Calculating and adding the area based on above information
        NUM_POINTS_THRESHOLD = args.noPointscluster
        FLAT_SC = []
        AREA_A2 = []
        AREA_P = []
        SEQ_COMP = []
        kdHydrophobicity = [] ## Kyte and Doolittle
        centroids = []
        i = 0

        for clus in range(no_binding_sites-1):
            seq_composition = []
            kd = 0.00
            cluster = torch.where(cluster_for_all_points==i+1)
            i+=1
            if len(cluster[0]) > NUM_POINTS_THRESHOLD:
                cluster2residue = point2residue[cluster[0]]
                INDEX = torch.unique(cluster2residue)

                clustercoords = residue2coords[INDEX]
                # calculate the centorid
                centroid_of_cluster = np.mean(clustercoords.numpy(), axis=0)                
                tmp = torch.tensor([[centroid_of_cluster[0], centroid_of_cluster[1], centroid_of_cluster[2]]])
                dist = torch.cdist(clustercoords, tmp, p=2)
                # calculate the closest point to the centroid
                centroid_point = clustercoords[torch.argmin(dist)]
                # calculate the closest residue
                distances_allvsall = torch.norm(centroid_point - residue2coords, dim=-1)
                distances_min = torch.argmin(distances_allvsall, dim=0)
                centroids.append(residue_number[distances_min.item()])

                area = 0
                for idx in INDEX:
                    area += sasa[idx+1]
                    seq_composition.append(f'{residue_number[idx]}_{sequence[idx]}')
                    ## Considering SASA < 20 as buried and >20 as exposed
                    ## To calculate the KD (hydrophobicity) we can sum up kd values
                    ## for the exposed residues.
                    if sasa[idx+1] > 20:
                        kd += float(AAKD2NUM[sequence[idx]])

                # print("area by summing up SASA:", area, "A^2") 
                # print("area by considering number of points:", len(cluster[0]))

                AREA_A2.append(area)
                AREA_P.append(len(cluster[0]))
                SEQ_COMP.append({"cluster_"+str(i+1): seq_composition})
                kdHydrophobicity.append(kd) ## save cluster number: [centeroid residue, kd value of the binding site], or just simply the kd
                FLAT_SC.append(flat_score(coords[cluster[0]])) ## return three numbers to be used to evaluate flatness.

        #print(acc, '\ncentroids: ', centroids,
        #        '\nflat_score: ', FLAT_SC,
        #        '\nKyte_Doolittle:', kdHydrophobicity)

        metrics.append({'acc': acc,
                        'sequence': sequence,
                        'sasa': sasa,
                        'point2residue': point2residue.numpy(),
                        'no_binding_sites': [no_binding_sites, len(centroids)],
                        'filtered_site_points': filtered_site_points.numpy(),
                        'normalized_labels': normalized_labels,
                        'cluster_for_all_points': cluster_for_all_points.numpy(),
                        'area_angstrom2': AREA_A2,
                        'area_points': AREA_P,
                        'sequence_compositions_all_clusters': SEQ_COMP,
                        'centroids': centroids,
                        'flat_score': FLAT_SC,
                        'Kyte_Doolittle': kdHydrophobicity})

        tt = cluster_for_all_points.numpy()
        norm_label_for_presentation = (tt-min(tt))/(max(tt)-min(tt))

        # make an empty file and then write it
        # with open(os.path.join(args.outpath,acc+"_pointcloud_clusters.pdb"), mode='w'): pass
        save_pointcloud(target_pdb, acc, coords, norm_label_for_presentation, args.outpath)

    pd.DataFrame(metrics).to_csv(os.path.join(args.data, 'result_metrics_117.csv'), index=False)
    # torch.save(metrics, os.path.join(args.data, 'result_metrics_all.pt')) # there is a error in storing sasa double vector.


if __name__ == '__main__':
    main(args=parser.parse_args())
