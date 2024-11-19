# This script is to analyse the output from soluble MPNN to select best candidates for AF2 filteration

## 1.0 Libraries
import os, glob, datetime, argparse, matplotlib
import subprocess
import pandas as pd

from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt


## 2.0 Functions
def create_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parse = argparse.ArgumentParser()
    
    parse.add_argument("--mpnn_folder" ,   type=str,   nargs=1,  required=True,  help="Path to folder containing MPNN run outputs.")
    parse.add_argument("--num_of_seqs" ,   type=int,   nargs=1,  required=True,  help="Number of MPNN sequences to fold.")
    parse.add_argument("--path_to_conda" , type=str,   nargs=1,  required=True,  help="Path to source conda.sh.")
    parse.add_argument("--path_to_env" ,   type=str,   nargs=1,  required=True,  help="Path to activate conda environment.")

    return parse

def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args

def plot_scores(df_parsed:pd.DataFrame, df_sel:pd.DataFrame, output:str):
    """  """
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[12,5])

    sns.violinplot(data=df_parsed, y='global_score', color='#DCDCDC', inner=None, ax=ax[0])
    sns.violinplot(data=df_sel,    y='global_score', color='#00BFFF', inner=None, ax=ax[0])
    sns.violinplot(data=df_parsed, y='score', color='#DCDCDC', inner=None, ax=ax[1])
    sns.violinplot(data=df_sel,    y='score', color='#32CD32', inner=None, ax=ax[1])
    sns.violinplot(data=df_sel,    y='len_seq', color='#FF8C00', inner=None, ax=ax[2])
    
    # Format labels
    ax[0].set_ylabel('pMPNN Global Score', fontsize=12, labelpad=10.0)
    ax[1].set_ylabel('pMPNN Score', fontsize=12, labelpad=10.0)
    ax[2].set_ylabel("Sel designs' lengths", fontsize=12, labelpad=10.0)
    
    plt.suptitle(f'{df_sel.shape[0]}/{df_parsed.shape[0]} Selected Sequences', fontsize=14)
    plt.tight_layout(h_pad=5., w_pad=5.)
    
    plt.savefig(f"{output}/{datetime.date.today()}_pmpnn_scores.svg", format='svg', dpi=200)
    plt.savefig(f"{output}/{datetime.date.today()}_pmpnn_scores.png", format='png', dpi=200)
    plt.close()
    return None

def main():
    """ Main Execution point """
    # Parse arguments
    args = parse_args(create_parser())
    root_folder = args.mpnn_folder[0]
    path_to_conda = args.path_to_conda[0]
    path_to_env = args.path_to_env[0]
    num_of_seqs = args.num_of_seqs[0]
    
    results_folder = [f.strip() for f in glob.iglob(os.path.join(root_folder, '*')) if "output" in f
                                                                                    if os.path.isdir(f)][0]
    optimisation_strategy = os.path.basename(results_folder)
    
    # make output dir
    output = f"{root_folder}/{optimisation_strategy}_analysis/"
    os.makedirs(output)

    # create folder for the fa files to submit to af2
    af2_monomer_fa = f'{output}/af2_monomer_fa_files/'
    os.makedirs(af2_monomer_fa)
    
    # create af2 pred output folder
    af2_monomer_pred = f'{output}/af2_monomer_pred/'
    os.makedirs(af2_monomer_pred)
    
    # create a folder for the slurm files
    slurm_output_folder = f'{output}/af_monomer_slurms/'
    os.makedirs(slurm_output_folder)

    # parse paths to fasta files
    des_paths = [f for f in glob.iglob(os.path.join(results_folder, '*/seqs/*.fa'))]
    print(f"{len(des_paths)} designs parsed!")

    # parse the fasta files
    #container
    parse_dic = {
        'design': [],
        'sequence':[], 
        'score':[],
        'global_score':[],
        'sequence_recovery':[]
        }

    for des_path in tqdm(des_paths, desc='Parsing Sol MPNN output'):
        header   = [line.strip() for line in open(des_path) if line.startswith('>')]
        sequence = [line.strip() for line in open(des_path) if not line.startswith('>')]

        for i in range(0, len(header), 1):
            des_name = [h.strip() for h in header[0].split(',')][0]
            if i == 0:
                parse_dic['design'].append(f'{des_name.replace("_","x")}_original')
                parse_dic['sequence'].append(sequence[i])
                parse_dic['score'].append(float([h.strip() for h in header[i].split(',')][1].split('=')[-1]))
                parse_dic['global_score'].append(float([h.strip() for h in header[i].split(',')][2].split('=')[-1]))
                parse_dic['sequence_recovery'].append(float(1))
            else:
                parse_dic['design'].append(f'{des_name.replace("_","x")}_mpnn_{[h.strip() for h in header[i].split(",")][1].split("=")[-1]}')
                parse_dic['sequence'].append(sequence[i])
                parse_dic['score'].append(float([h.strip() for h in header[i].split(',')][2].split('=')[-1]))
                parse_dic['global_score'].append(float([h.strip() for h in header[i].split(',')][3].split('=')[-1]))
                parse_dic['sequence_recovery'].append(float([h.strip() for h in header[i].split(',')][4].split('=')[-1]))
    df_parsed = pd.DataFrame(parse_dic)
    df_parsed.to_csv(f'{output}/all_sol_mpnn_designs.csv')

    # drop duplicated sequences
    df_uniq = df_parsed.copy().sort_values(by='global_score').drop_duplicates(subset='sequence', keep='first').reset_index(drop=True)
    df_uniq['len_seq'] = df_uniq['sequence'].apply(lambda x: len(x))
    df_uniq.sort_values(by='len_seq',ascending= True).reset_index(drop=True)

    # Get top 5 mpnn designs for further analysis based on lowest global score
    #container
    vessel = []

    #add the design number to group with
    df_uniq['parent'] = df_uniq['design'].apply(lambda x: x.split('_')[0].lstrip('>'))

    for g in tqdm(df_uniq.groupby(by='parent')):
        df_temp = g[1]
        df_temp.sort_values(by='global_score', ascending=True, inplace=True)
        df_temp.reset_index(drop=True, inplace=True)
        vessel.append(df_temp.head(num_of_seqs))
    df_sel = pd.concat(vessel).reset_index(drop=True)
    df_sel.to_csv(f'{output}/selected_sol_mpnn_designs.csv')

    # plot the scores
    plot_scores(df_parsed=df_parsed, df_sel=df_sel, output=output)

    # write out the fasta files for the AF2 prediction
    # sort the files by length (repeated just to make sure)
    df_sel.sort_values(by='len_seq', ascending=True, inplace=True)
    df_sel.reset_index(drop=True, inplace=True)
    
    # split the list to 5 GPUs with extended time
    interval = int(df_sel.shape[0]/5)
    
    # write out the fasta files
    for start in tqdm(range(0, (df_sel.shape[0]), interval), desc='Writing FASTA files'):
        with open(f'{af2_monomer_fa}/af_monomer_{start}.fa', 'w') as w:
            for i,label,seq in zip(df_sel.index, df_sel.design, df_sel.sequence):
                if i >= start and i <= start+interval:
                    w.write(f'{label}' + '\n')
                    w.write(f'{seq}' + '\n')

    # export the list of paths to the fasta files
    fasta_paths = [os.path.abspath(x) for x in glob.iglob(os.path.join(af2_monomer_fa, '*.fa'))]
    with open(f'{output}/fasta_paths.list', 'w') as w:
        for fp in range(0, len(fasta_paths), 1):
            w.write(fasta_paths[fp] + '\n')
            
            
    # write the slurm files
    for fasta_file in tqdm(fasta_paths, desc='Writing slurm files'):
        # set requested CPU time
        slurm_file = """#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem 80gb
#SBATCH --time 12:00:00
#SBATCH --output=colabfold_monomer_%A.log
# load the correct environment 
source {path2conda}
conda activate {path2env}
SECONDS=0
# load modules
module purge
module load gcc/11.3.0
module load cuda/11.8.0
module load cudnn/8.7.0.84-11.8
#Run prediction
INPUT={input_fa}
OUTPUT={output_folder}
colabfold_batch $INPUT \
$OUTPUT \
--msa-mode single_sequence \
--num-recycle 3 \

printf '%s\n' '------------------------------------------------'
ELAPSED="Prediction took: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
printf $ELAPSED
        """.format(
            path2conda=path_to_conda,
            path2env=path_to_env,
            input_fa=fasta_file,
            output_folder=af2_monomer_pred)
        
        # make the slurm file
        with open(f"{slurm_output_folder}/{os.path.basename(fasta_file).rstrip('.fa')}.slurm", 'w') as w01:
            w01.write(slurm_file)
            w01.close()
            
    # submit slurm files
    for slurm_ in glob.iglob(os.path.join(slurm_output_folder, '*.slurm')):
        subprocess.call(['sbatch', slurm_])
        
    # Clean up
    try:
        [f.strip() for f in glob.iglob(os.path.join(root_folder, '*')) if "input" in f if os.path.isdir(f)][0]
    except IndexError:
        print('Sol MPNN input folder is compressed')
    else:
        sol_mpnn_folder_path = [f.strip() for f in glob.iglob(os.path.join(root_folder, '*')) if "input" in f if os.path.isdir(f)][0]
        print('Compressing Sol MPNN input folder')
        os.system(f'tar -czf {os.path.dirname(sol_mpnn_folder_path)}/{os.path.basename(sol_mpnn_folder_path)}.tar.gz {sol_mpnn_folder_path}')
        print("Removing uncompressed Sol MPNN input folder")
        os.system(f"rm -r {sol_mpnn_folder_path}")
    
    return None

# Execute
if __name__ == "__main__":
    main()
