# This script is to analyse the the af2 prediction for MPNN sequences
# this is adapted to work with colabfold 1.5.0

## 1.0 Libraries
import os, glob, shutil, argparse, matplotlib
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm

# import plotting
import seaborn as sns
import matplotlib.pyplot as plt

# import Biopython
from Bio.PDB import *
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis


## 2.0 Functions
def create_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parse = argparse.ArgumentParser()
    
    parse.add_argument("--mpnn_folder",            type=str,   nargs=1,  required=True,  help="Path to ProteinMPNN run output folder.")
    parse.add_argument("--rs_models_dir",          type=str,   nargs=1,  required=True,  help="Path to folder containing rosetta designed models.")
    parse.add_argument("--rs_model_binder_chain",  type=str,   nargs=1,  required=True,  help="Rosetta model binder chain ID.")
    parse.add_argument("--af_model_binder_chain",  type=str,   nargs=1,  required=True,  help="AlphaFold2 model chain ID.")
    parse.add_argument("--pLDDT_thresh"         ,  type=float, nargs=1,  required=False, help="pLDDT threshold for filtering, Default 80.")
    parse.add_argument("--RMSD_thresh"          ,  type=float, nargs=1,  required=False, help="C⍺ RMSD to model threshold for filtering, Default 1.5Å.")
    return parse

def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args

def calculate_rmsd(rs_model:str, af_model:str, rmsd_out:str, rs_model_binder_chain:str='B', af_model_binder_chain:str='A'):
    
    # parse the rs model
    designname = os.path.basename(rs_model).replace('.pdb','')
    parser1 = PDBParser(QUIET=True)
    struct = parser1.get_structure(designname, rs_model)
    res_binder_rs = Selection.unfold_entities(struct, 'R')
    ca_atoms_binder_rs = [x['CA'] for x in res_binder_rs if x.get_parent().get_id() in rs_model_binder_chain]
    
    # parse the af modes
    af_model_name = os.path.basename(af_model).replace('.pdb','')
    parser2 = PDBParser(QUIET=True)
    structaf2 = parser2.get_structure(af_model_name, af_model)
    moving_atoms = Selection.unfold_entities(structaf2, 'A')
    res_target_af2= Selection.unfold_entities(structaf2, 'R')
    ca_atoms_binder_af2 = [x['CA'] for x in res_target_af2 if x.get_parent().get_id()== af_model_binder_chain] 
    
    # align the 2 models based on CA of the target
    superimposer = Superimposer()
    superimposer.set_atoms(ca_atoms_binder_rs, ca_atoms_binder_af2)
    superimposer.apply(moving_atoms)
    
    # OUTPUT  
    # af model
    pdbio = PDBIO()
    pdbio.set_structure(structaf2)
    pdbio.save(f'{rmsd_out}/{os.path.basename(af_model).replace(".pdb","")}.pdb')
    
    # compute the rmsd for binder
    design_binder_ca_coord = np.array([x.get_coord() for x in ca_atoms_binder_rs])
    af2_binder_ca_coord = np.array([x.get_coord() for x in ca_atoms_binder_af2])  
    
    # calculate rmsd of binder
    rmsd = design_binder_ca_coord - af2_binder_ca_coord
    rmsd = np.square(rmsd)
    rmsd = np.sum(rmsd, axis=1)
    rmsd = np.mean(rmsd)
    rmsd = np.sqrt(rmsd)    
    
    return (designname, rmsd, af_model_name)

def binder_properties_calculator_from_df(df_prop:pd.DataFrame, sequence_column:str):
    """
    This function will use the sequences in the df
    ExPASY protparam to caluclate the grafts' properties
    
    - df_prop: the DataFrame that will be used for the analysis
    - item_for_sorting: the column name that will be used to sort the values in the resulted DF
    - sequence_column: the name of the column containing the sequnce information for the designs
    - outpath: the dir where the resulted df will be exported
    - identifier: the desired name to identify the run
    """
    
    df_prop["Mwt"]  = df_prop[f'{sequence_column}'].apply(lambda x: round(ProteinAnalysis(x).molecular_weight(),2))
    df_prop["IsoelectricPoint"] = df_prop[f'{sequence_column}'].apply(lambda x: round(ProteinAnalysis(x).isoelectric_point(),2))
    df_prop["charge_at_pH_7.4"] = df_prop[f'{sequence_column}'].apply(lambda x: round(ProteinAnalysis(x).charge_at_pH(7.4),2))
    df_prop["GRAVY"] = df_prop[f'{sequence_column}'].apply(lambda x: round(ProteinAnalysis(x).gravy(),2))
    df_prop["instability_index"] = df_prop[f'{sequence_column}'].apply(lambda x: round(ProteinAnalysis(x).instability_index(),2))
    df_prop["aromaticity"] = df_prop[f'{sequence_column}'].apply(lambda x: round(ProteinAnalysis(x).aromaticity(),2))
        
    #reset index
    df_prop.reset_index(drop=True, inplace=True)
           
    return df_prop

def plot_prop(df_filtered:pd.DataFrame, df_sel:pd.DataFrame, outpath:str, iden:str, figsize:list, export:bool=True):
    """
    """
    fig, axs = plt.subplots( nrows=1, ncols=6, figsize=figsize )

    sns.violinplot(y='binder_rmsd',      data=df_filtered, ax=axs[0], color='#00FF7F')
    sns.violinplot(y='score',            data=df_filtered, ax=axs[1], color='#00BFFF')
    sns.violinplot(y='global_score',     data=df_filtered, ax=axs[2], color='#F4A460')
    sns.violinplot(y='GRAVY',            data=df_filtered, ax=axs[3], color='#7FFFD4')
    sns.violinplot(y='charge_at_pH_7.4', data=df_filtered, ax=axs[4], color='#DA70D6')
    sns.violinplot(y='instability_index',data=df_filtered, ax=axs[5], color='#DCDCDC', inner=None)
    sns.violinplot(y='instability_index',data=df_sel,      ax=axs[5], color='#9370DB')

    axs[0].set_ylabel('Binder RMSD (Å)', fontsize=12, labelpad=5.)
    axs[1].set_ylabel('MPNN score', fontsize=12, labelpad=5.)
    axs[2].set_ylabel('MPNN global score', fontsize=12, labelpad=5.)
    axs[3].set_ylabel('GRAVY index', fontsize=12, labelpad=5.)
    axs[4].set_ylabel('Charge at pH 7.4', fontsize=12, labelpad=5.)
    axs[5].set_ylabel('Instability Index', fontsize=12, labelpad=5.)

    plt.suptitle(f"""{df_sel.shape[0]} selected designs
                 """, fontsize=22)
    plt.tight_layout(h_pad=5.0, w_pad=3.0)

    if export:
        plt.savefig(f'{outpath}/{iden}.svg', format="svg" ,dpi=300, transparent=True)
    
    plt.close()
    return None

def plot_all_decoys(df_sel01:pd.DataFrame, df_avg:pd.DataFrame, output:str):
    """  """
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[8,6])

    sns.violinplot(data=df_avg,    y='avg_ptm',               ax=axs[0],    color='#DCDCDC', inner=None)
    sns.violinplot(data=df_sel01,  y='avg_ptm',               ax=axs[0],    color='#00BFFF', inner=None)
    sns.violinplot(data=df_avg,    y='avg_binder_mean_plddt', ax=axs[1],    color='#DCDCDC', inner=None)
    sns.violinplot(data=df_sel01,  y='avg_binder_mean_plddt', ax=axs[1],    color='#00FF7F', inner=None)

    # format labels
    axs[0].set_ylabel('PTM score (averaged)',   fontsize=12, labelpad=10.)
    axs[1].set_ylabel('pLDDT score (averaged)', fontsize=12, labelpad=10.)

    plt.suptitle(f"{df_sel01.shape[0]} / {df_avg.shape[0]} selected decoys", fontsize=18, )
    plt.tight_layout(w_pad=5.0, h_pad=10.0)

    plt.savefig(f'{output}/all_decoys_params.svg', format='svg', dpi=200, transparent=True)
    plt.close()
    return None

def af2_json_parser(json_paths:list, output:str) -> pd.DataFrame:
    """ This function is to parse af2 json files into a dataframe """
    # container
    info_dic = {'design':[], 'ptm':[], 'max_pae':[], 'binder_plddt':[]}
    # get the info
    for json_file in tqdm(json_paths, desc='Parsing AF2 scores'):
        df_temp = pd.read_json(json_file)
        info_dic['design'].append(os.path.basename(json_file).replace('.json','').replace('scores','unrelaxed'))
        info_dic['ptm'].append(float(df_temp.ptm.to_list()[0]))
        info_dic['max_pae'].append(float(df_temp.max_pae.to_list()[0]))
        info_dic['binder_plddt'].append(df_temp.plddt.to_list())
    
    df_parsed = pd.DataFrame(info_dic)
    df_parsed['binder_mean_plddt'] = df_parsed['binder_plddt'].apply(lambda x: np.mean(x))
    df_parsed['src'] = df_parsed['design'].apply(lambda x: x.split('_unrelaxed')[0])
    df_parsed.sort_values(by='ptm', ascending=False, inplace=True) # sort values by ptm score
    df_parsed.reset_index(drop=True, inplace=True) 
    df_parsed.to_csv(f'{output}/af2_monomer_all_decoys_analyzed.csv')
    return df_parsed

def avg_metrics(df_parsed:pd.DataFrame, output:str) -> pd.DataFrame:
    """ This function is to average the metrics from all 5 AF2 models for designs filteration"""

    # container
    vessel01 = {'label':[], 'design':[], 'avg_ptm':[], 'avg_max_pae':[], 'avg_binder_mean_plddt' :[]}
    
    for gp01 in tqdm(df_parsed.groupby(by='src'), desc='Averaging models scores'):
        df_temp01 = gp01[1].reset_index(drop=1)
        df_temp01.sort_values(by='binder_mean_plddt', ascending=False, inplace=True)
        df_temp01.reset_index(drop=True, inplace=True)
        vessel01['label'].append(f">{gp01[0]}")
        vessel01['design'].append(df_temp01.design[0])
        vessel01['avg_ptm'].append(df_temp01.ptm.mean())
        vessel01['avg_max_pae'].append(df_temp01.max_pae.mean())
        vessel01['avg_binder_mean_plddt'].append(df_temp01.binder_mean_plddt.mean())
    df_avg = pd.DataFrame(vessel01)
    df_avg.to_csv(f'{output}/af2_monomer_all_analyzed_decoys_averaged.csv')
    return df_avg

def extract_target_info(path:str, target_pdb_id:str, target_chain:str, output:str):
    """ This function is to extract the target protein """
    # load target
    parser1 = PDBParser(QUIET=True)
    struct = parser1.get_structure(os.path.basename(path), path)[0][target_chain]
    
    # export target
    pdbio = PDBIO()
    pdbio.set_structure(struct)
    pdbio.save(f'{output}/{target_pdb_id}.pdb')
    
    print(f"{target_pdb_id} template exported !")
    
    target_seq = SeqUtils.seq1("".join([r.get_resname() for r in struct.get_residues()]))
    return target_seq

def main():
    """ Main Execution point """

    # Parse arguments
    args = parse_args(create_parser())
    root_folder = args.mpnn_folder[0]
    rs_models_folder = args.rs_models_dir[0]
    rs_binder_chain = args.rs_model_binder_chain[0]
    af2_binder_chain = args.af_model_binder_chain[0]

    if args.pLDDT_thresh != None:
        pLDDT_thresh = float(args.pLDDT_thresh[0])
    else:
        pLDDT_thresh = 80.0

    if args.RMSD_thresh != None:
        RMSD_thresh = float(args.RMSD_thresh[0])
    else:
        RMSD_thresh = 1.5

    #### Get the Needed Paths ####

    ##### path to the rosetta models
    rs_models_paths  = [f.strip() for f in glob.iglob(os.path.join(rs_models_folder, '*.pdb'))]

    ##### path to the mpnn output analysis 
    mpnn_output_analysis_dir = [f.strip() for f in glob.iglob(os.path.join(root_folder, '*')) if '_output_analysis' in f
                                                                                              if os.path.isdir(f) ][0]
    ##### path to the mpnn output score file
    mpnn_sc = [f.strip() for f in glob.iglob(os.path.join(mpnn_output_analysis_dir,'*.csv')) if 'selected' in f][0]

    ##### path to the af2 monomer predictions 
    af2_monomer_pred_dir = [f.strip() for f in glob.iglob(os.path.join(mpnn_output_analysis_dir, '*')) if 'monomer_pred' in f
                                                                                                       if os.path.isdir(f)][0]
    # parse AF2 prediction file paths
    json_paths  = [f for f in glob.iglob(os.path.join(af2_monomer_pred_dir, '*.json')) if 'scores' in f]

    #### Make output dirs ####
    # make output dir for analysis
    output = os.path.abspath(os.path.join(af2_monomer_pred_dir, "..", "af2_monomer_pred_analysis"))
    os.makedirs(output, exist_ok=0)

    # make output for the RMSD files
    rmsd_out = os.path.join(output, "aligned_models")
    os.makedirs(rmsd_out, exist_ok=0)

    # make dir for selected designs
    sel_des_out = os.path.join(output, "sel_designs")
    os.makedirs(sel_des_out, exist_ok=0)

    ### Parse AF output and select best decoys ###
    df_parsed = af2_json_parser(json_paths=json_paths, output=output)
    print(f"{df_parsed.shape[0]} parsed decoys")

    # Add the averages for all models to seek AF2 agreement within models
    df_avg = avg_metrics(df_parsed=df_parsed, output=output)
    print(f'{df_avg.shape[0]} averaged decoys')

    # filter binders based on ptm and binder plddt
    df_sel01 = df_avg[(df_avg.avg_binder_mean_plddt >= pLDDT_thresh)]
    df_sel01.sort_values(by='avg_binder_mean_plddt', ascending=False, inplace=True)
    df_sel01.reset_index(drop=1, inplace=True)
    df_sel01.to_csv(f'{output}/af2_monomer_best_decoys.csv')
    assert df_sel01.shape[0] != 0, "No decoys passed the filters"
    print('af2_monomer_best_decoys Exported!')

    # plt all decoys
    plot_all_decoys(df_sel01=df_sel01, df_avg=df_avg, output=output)

    #get the paths to the sel pdbs
    sel_des_names = [d for d in df_sel01.design]
    sel_des_paths = [f'{af2_monomer_pred_dir}/{d}.pdb' for d in sel_des_names]
    
    # map the rosetta models to their respective af2 prediction of the optimised model
    models_map = {}
    for rs in tqdm(rs_models_paths, desc='Mapping AF2 models to Rs models'):
        rs_name = os.path.basename(rs).replace('.pdb','')
        af_correspondence = [a for a in sel_des_paths if rs_name == os.path.basename(a).split('_mpnn')[0].replace("x","_")]
        if len(af_correspondence) > 0:
            models_map[rs] = af_correspondence
    assert len(models_map.keys()) != 0, "Model mapping failure"

    # calculate the rmsd between the rs models and the af models
    rmsd_dic = {'design':[], 'binder_rmsd':[], 'af2_model':[]}

    for key, value in tqdm(models_map.items(), desc='Calculating RMSD'):
        rs_model_path = key
        for v in value:
            af_model_path = v
            try:
                designname, rmsd, af_model_name = calculate_rmsd(rs_model=rs_model_path, af_model=af_model_path, rmsd_out=rmsd_out, rs_model_binder_chain=rs_binder_chain, af_model_binder_chain=af2_binder_chain)
            except:
                print(f'{key}, {v} not matched!!')
            else:
                rmsd_dic['design'].append(designname)
                rmsd_dic['binder_rmsd'].append(rmsd)
                rmsd_dic['af2_model'].append(af_model_name)
    df_rmsd = pd.DataFrame(rmsd_dic)
    df_rmsd.to_csv(os.path.join(output, "af2_monomer_folded_decoys_rmsd.csv"))

    # select promising sequences that folded correctly with CA RMSD <= 1.5
    df_sel02 = df_rmsd[(df_rmsd.binder_rmsd <= RMSD_thresh)].sort_values(by='binder_rmsd', ascending=1).reset_index(drop=True)
    df_sel02['label'] = df_sel02['af2_model'].apply(lambda x: f'>{x.split("_unrelaxed")[0]}')
    df_sel02.drop(columns=['design'], inplace=True)
    df_sel02.to_csv(os.path.join(output,"af2_monomer_filtered_decoys_by_rmsd.csv"))

    # Add the binder RMSD information to the selected designs
    df_merged01 = pd.merge(left=df_sel01, right=df_sel02, on='label', how='inner')
    df_merged01.reset_index(drop=True, inplace=True)

    # load the MPNN info 
    df_info = pd.read_csv(mpnn_sc, index_col=0).rename(columns={'design':'label'})
    df_merged02 = pd.merge(left=df_info, right=df_merged01, on='label', how='inner')

    # calculate the sequence properties
    df_prop = binder_properties_calculator_from_df(df_prop=df_merged02, sequence_column="sequence")
    df_prop['design_strategy'] = os.path.basename(root_folder).split('mpnn_')[-1]
    df_prop.to_csv(os.path.join(output, "af_monomer_selection.csv"))

    # seperate the selected designs 
    for model in tqdm(df_merged02.af2_model, desc='Copying selected AF2 models'):
        shutil.copy(f'{rmsd_out}/{model}.pdb', f'{sel_des_out}/{model}.pdb')

        
    return None

# Execute
if __name__ == "__main__":
    main()
