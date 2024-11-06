#!/usr/bin/env python
"""
Created on Thu Feb  16 11:28:23 2023

Modified on Fri Feb 17 17:11:00 2023

@author: hamedkhakzad
"""

## clustering masif seed results

import argparse
import matplotlib
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
import numpy as np
import seaborn as sns
import pandas as pd
import os
from Bio.PDB import *
import logomaker


def read_inp(site):

    base_dir = 'analysis_fixed_size_site_{}/'.format(site)
    data_dir = os.path.join(base_dir, 'out_data/')
    # Read the list of results.
    with open(os.path.join(base_dir, 'list_out_peptides_site_{}.txt'.format(site))) as f:
        list_names = f.readlines()
        list_names = [x.rstrip() for x in list_names]

    # Read all of the pairwise rmsd distances 
    mydata_dict = {}

    for i in range(len(list_names)):
        fn = 'rmsd_{}.npy'.format(i)
        try:
            mydata_dict[int(i)] = np.load(os.path.join(data_dir,fn))
        except:
            continue
    return base_dir, mydata_dict, list_names

def remove_error(mydata_dict, list_names_1):

    ## Remove any values that had an error.
    active_keys = list(mydata_dict.keys())
    mydata = []
    for key in active_keys:
        vals = np.array(mydata_dict[key])
        vals = vals[active_keys]
        mydata.append(vals)
    mydata = np.asarray(mydata)
    # remove those results that had problems.
    list_names = [x for ix, x in enumerate(list_names_1) if ix in active_keys]
    
    return mydata, list_names

def plot_rmsd(base_dir, mydata):

    # Plot the rmds. 
    mydata_flat = np.reshape(mydata, [-1])
    plt.hist(mydata_flat)
    plt.savefig(os.path.join(base_dir, 'rmsd.pdf'), type='pdf')
    
    return

def mds_clustering(base_dir, mydata):

    # Run MDS on all pairwise distances, 
    runmds = MDS(n_components=2, metric=True, dissimilarity="precomputed")
    mds_plot = runmds.fit_transform(mydata)
    np.save(os.path.join(base_dir, 'mds_plot.npy'), mds_plot)

    # Rotate the mds_plot by 180 degrees
    rot_matrix = [[-1, 0], [0,-1]]
    mds_plot = np.dot(mds_plot, rot_matrix)
    
    kmeans = KMeans(n_clusters=6, random_state=0).fit(mds_plot)
    myplot = pd.DataFrame(mds_plot, columns=['x', 'y'])
    myplot['cluster'] = kmeans.labels_

    # plot the clsuters
    histpy = np.histogram(kmeans.labels_, bins=[0,1,2,3,4,5,6])
    histpy = np.array([histpy[0], [0,1,2,3,4,5]])

    histdata = pd.DataFrame(histpy.T, columns=['Count', 'Cluster number'])
    sns.barplot(data=histdata,  y='Count', x='Cluster number')
    plt.savefig(os.path.join(base_dir, 'barplot.pdf'), type='pdf')

    return mds_plot, kmeans.labels_, kmeans.cluster_centers_

def plot_cluster(base_dir, mds_plot, labels, centers):

    palette=sns.color_palette('muted', n_colors=6)

    myplot = pd.DataFrame(mds_plot, columns=['x', 'y'])
    myplot['cluster'] = labels
    plt.figure(figsize=(10,8))
    ax = sns.scatterplot(data=myplot, x="x", y='y', hue='cluster', palette=palette)
    markers = {"star": "*", 'triangle': '>', 'hexagon': 'h'}
    ax.set_xlabel('')
    ax.set_ylabel('')
    # Mark the patch of structures selected.
    cc = centers
    for k in range(len(cc)):
        mypatch1 = matplotlib.patches.Rectangle((cc[k][0]-0.75,cc[k][1]-0.75), width=1.5, height=1.5, color='gray', fill=True, alpha=0.25)
        mypatch2 = matplotlib.patches.Rectangle((cc[k][0]-0.75,cc[k][1]-0.75), width=1.5, height=1.5, color='black', fill=False)
        ax.add_patch(mypatch1)
        ax.add_patch(mypatch2)

    plt.savefig(os.path.join(base_dir, 'cluster_by_rmsd_legend.pdf'), format='pdf')

    return

def plot_logo(base_dir, list_names, mds_plot, centers):

    ## Draw a logoplot on the selected helices positions of the cluster of choice.
    # Top cluster!
    top_cluster = 1
    # Find points near cluster center top_cluster
    cp = centers[top_cluster]
    idx_cluster_1_box = []
    for i in range(len(mds_plot)):
        p = mds_plot[i]
        dist1 = np.sqrt(np.sum(np.square(p - cp)))
        if dist1 < 0.75:
            idx_cluster_1_box.append(i)

    parser = PDBParser()

    # One hot encoding of the residues.
    data_one_hot = np.zeros((12, 20))
    col_names = [Polypeptide.index_to_one(x) for x in range(20)]
    pdb_dir = os.path.join(base_dir, 'out_pdb')

    for i in range(len(list_names)):
        if i in idx_cluster_1_box:
            compstruct = parser.get_structure('comp', os.path.join(pdb_dir,list_names[i]+'_frag.pdb'))
            compres = Selection.unfold_entities(compstruct, 'R')
            do_compare=False
            
            rrix = 0        

            resn0 = Polypeptide.three_to_one(compres[rrix].get_resname())
            resn1 = Polypeptide.three_to_one(compres[rrix+1].get_resname())
            resn2 = Polypeptide.three_to_one(compres[rrix+2].get_resname())
            resn3 = Polypeptide.three_to_one(compres[rrix+3].get_resname())
            resn4 = Polypeptide.three_to_one(compres[rrix+4].get_resname())
            resn5 = Polypeptide.three_to_one(compres[rrix+5].get_resname())
            resn6 = Polypeptide.three_to_one(compres[rrix+6].get_resname())
            resn7 = Polypeptide.three_to_one(compres[rrix+7].get_resname())
            resn8 = Polypeptide.three_to_one(compres[rrix+8].get_resname())
            resn9 = Polypeptide.three_to_one(compres[rrix+9].get_resname())
            resn10 = Polypeptide.three_to_one(compres[rrix+10].get_resname())
            resn11 = Polypeptide.three_to_one(compres[rrix+11].get_resname())
            
            # Store the residue in the logo plot.
            data_one_hot[0][col_names.index(resn0)] += 1
            data_one_hot[1][col_names.index(resn1)] += 1
            data_one_hot[2][col_names.index(resn2)] += 1
            data_one_hot[3][col_names.index(resn3)] += 1
            data_one_hot[4][col_names.index(resn4)] += 1
            data_one_hot[5][col_names.index(resn5)] += 1
            data_one_hot[6][col_names.index(resn6)] += 1
            data_one_hot[7][col_names.index(resn7)] += 1
            data_one_hot[8][col_names.index(resn8)] += 1
            data_one_hot[9][col_names.index(resn9)] += 1
            data_one_hot[10][col_names.index(resn10)] += 1
            data_one_hot[11][col_names.index(resn11)] += 1    
    
    data_one_hot_norm = data_one_hot/np.sum(data_one_hot[0])
    df_aligns = pd.DataFrame(data_one_hot_norm, columns = col_names)

    ## Making logo
    fig, ax = plt.subplots(1,1,figsize=[20,10])

    # create Logo object
    crp_logo = logomaker.Logo(df_aligns, ax=ax, color_scheme='skylign_protein', width=0.8)
    # style using Logo methods
    #crp_logo.style_spines(visible=False)
    crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
    crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    #crp_logo.ax.set_ylabel("Frequency")
    crp_logo.ax.xaxis.set_ticks_position('none')
    crp_logo.ax.xaxis.set_tick_params(pad=-1)

    plt.savefig(os.path.join(base_dir, 'logoplot.pdf'))

    ax = sns.heatmap(df_aligns.transpose())
    ax.set_yticklabels(col_names, rotation=0)
    plt.savefig(os.path.join(base_dir, 'heatmap_seed.pdf'))

    return


def main():
        
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--site", help="the site number to target", default=0, type=int)
    args = parser.parse_args()
    
    site = args.site

    base_dir, mydata_dict, list_names_1 = read_inp(site)
    mydata, list_names = remove_error(mydata_dict, list_names_1)
    plot_rmsd(base_dir, mydata)
    mds_plot, labels, centers = mds_clustering(base_dir, mydata)
    plot_cluster(base_dir, mds_plot, labels, centers)
    plot_logo(base_dir, list_names, mds_plot, centers)

if __name__ == '__main__':
    main()


