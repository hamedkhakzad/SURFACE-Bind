# this script is to put each design in a dedicated folder preparing them for processing with Soluble ProteinMPNN
# Also a list for the made folders' paths will be exported for usage ease

## 1.0 libraries
import os, glob, shutil, argparse
import subprocess
from tqdm import tqdm


## 2.0 Functions
def create_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parse = argparse.ArgumentParser()
    
    parse.add_argument("--rs_models_dir" ,         type=str,   nargs=1,  required=True,                                                   help="Path to folder containing rosetta designed models that will be input for optimisation.")
    parse.add_argument("--pmpnn_repo" ,            type=str,   nargs=1,  required=True,                                                   help="Path to ProteinMPNN repository.")
    parse.add_argument("--scripts_src_repo" ,      type=str,   nargs=1,  required=False,                                                  help="Path to scripts src folder.")
    parse.add_argument("--binder_chain" ,          type=str,   nargs=1,  required=True,                                                   help="Binder chain ID to design.")
    parse.add_argument("--interface_threshold" ,   type=float, nargs=1,  required=False,                                                  help="Distance threshold for a residue to be considered a hotspot,, Default 3.5A.")
    parse.add_argument("--seqs_to_generate" ,      type=int,   nargs=1,  required=True,                                                   help="Number of sequences to generate per design.")
    parse.add_argument("--path_to_conda" ,         type=str,   nargs=1,  required=True,                                                   help="Path to source conda.sh.")
    parse.add_argument("--path_to_env" ,           type=str,   nargs=1,  required=True,                                                   help="Path to activate conda environment.")
    parse.add_argument("--optimisation_strategy",  type=str,   nargs=1,  required=True,  choices=['fixed_interface', 'design_interface'], help="Name of ProteinMPNN strategy to be used, ['fixed_interface', 'design_interface'].")

    return parse

def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args

## 3.0 Execution
def main():
    "Main execution point"
    
    # Parse arguments
    args = parse_args(create_parser())
    pmpnn_repo = args.pmpnn_repo[0]
    binder_chain = args.binder_chain[0]
    rs_models_dir = args.rs_models_dir[0]
    seqs_to_generate = args.seqs_to_generate[0]
    path_to_conda = args.path_to_conda[0]
    path_to_env = args.path_to_env[0]
    optimisation_strategy = args.optimisation_strategy[0]
    
    if args.scripts_src_repo != None:
        scripts_src_repo = args.scripts_src_repo[0]
    
    if args.interface_threshold != None:
        interface_threshold = args.interface_threshold[0]
    
    # get paths for rosetta models to be analysed
    rs_models_paths = [f.strip() for f in glob.iglob(os.path.join(rs_models_dir, '*.pdb'))]
    
    # make output folder 
    output_folder = os.path.join(os.path.abspath("."), f"sol_mpnn_{optimisation_strategy}")
    os.makedirs(output_folder)
    
    # make sol mpnn input dir
    sol_mpnn_input_dir = os.path.join(output_folder, "sol_mpnn_input")
    os.makedirs(sol_mpnn_input_dir)
    
    # make sol mpnn output dir
    sol_mpnn_output_dir_root = os.path.join(output_folder, f"sol_mpnn_{optimisation_strategy}_output")
    os.makedirs(sol_mpnn_output_dir_root)
    
    # create a folder to export slurm files
    slurm_output_folder = os.path.join(output_folder, "slurms")
    os.makedirs(slurm_output_folder)
    
    # create output folder for temp_files
    _temp = os.path.join(output_folder, "_temp")
    os.makedirs(_temp)
    
    # prepare soluble mpnn input
    if optimisation_strategy == 'fixed_interface':
        for rs_model01 in tqdm(rs_models_paths, desc='Preparing Sol MPNN input, making folders'):
            new_folder01 = f"{sol_mpnn_input_dir}/{os.path.basename(rs_model01).replace('.pdb','')}/"
            os.makedirs(new_folder01)
            shutil.copy(rs_model01, os.path.join(new_folder01, os.path.basename(rs_model01)))   
    else:
        for idx01 in range(0,len(rs_models_paths), 200):
            new_folder02 = f"{sol_mpnn_input_dir}/slice_{idx01:04}"
            os.makedirs(new_folder02)
            for element01 in rs_models_paths[idx01:idx01+200]:
                shutil.copy(element01, os.path.join(new_folder02, os.path.basename(element01)))
            
    # get the design's path lists that will be passed to the pMPNN script
    pmpnn_input = [os.path.abspath(f.strip()) for f in glob.iglob(os.path.join(sol_mpnn_input_dir, '*'))]
    
    # export the full list for check up
    with open(f'{output_folder}/all_sol_mpnn_input_des_paths_{optimisation_strategy}.txt', 'w') as w01:
        for line in pmpnn_input:
            w01.write(line + '\n')
        w01.close()

    # write the slurm files
    for rs_model02 in tqdm(pmpnn_input, desc='writing slurm files'):
        # make an output dir for each design
        sol_mpnn_output_dir = f"{sol_mpnn_output_dir_root}/{os.path.basename(rs_model02).replace('.pdb','')}_out/"
        os.makedirs(sol_mpnn_output_dir)

        # set requested CPU time
        if optimisation_strategy == 'fixed_interface':
            request_time = '00:30:00'
        else:
            request_time = '8:00:00'
            
        slurm_head = """#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time {time}
#SBATCH --qos=serial

# load the correct environment 
source {path2conda}
conda activate {path2env}

# make a temp copy of the mpnn scripts
TEMP="{t}/{d}_temp/"
mkdir $TEMP
cp -r {r}/soluble_model_weights/ $TEMP/soluble_model_weights/
cp -r {r}/helper_scripts $TEMP/helper_scripts
cp {r}/protein_mpnn_run.py $TEMP/protein_mpnn_run.py
cp {r}/protein_mpnn_utils.py $TEMP/protein_mpnn_utils.py

#SET ROOT OF DIRECTORY
paramdir="$TEMP/soluble_model_weights/"
runpath="$TEMP/protein_mpnn_run.py"
helperscripts="$TEMP/helper_scripts"

path_for_parsed_chains={out}"/parsed_pdbs.jsonl"
path_for_assigned_chains={out}"/assigned_pdbs.jsonl"
path_for_fixed_positions={out}"/fixed_pdbs.jsonl"
        
        """.format(path2conda=path_to_conda,
            path2env=path_to_env,
            time=request_time,
            p=os.path.abspath(output_folder).split('/')[4],
            r=pmpnn_repo,
            t=_temp,
            d=os.path.basename(rs_model02).replace('.pdb',''),
            out=sol_mpnn_output_dir)
        
        if optimisation_strategy == 'fixed_interface':
            slurm_hotspots=""" 
#The first amino acid in the chain corresponds to 1 and not PDB residues index for now. 
HOTSPOTS=$(python3 {src}/find_hotspots.py {pdb_folder} {binder_chain} {threshold}) 
fixed_positions=$HOTSPOTS 
echo $fixed_positions 
        """.format(src=scripts_src_repo,
                      pdb_folder=rs_model02,
                      binder_chain=binder_chain,
                      threshold=interface_threshold)
        else:
            slurm_hotspots="fixed_positions=''"
        
        slurm_body = """
python3 $helperscripts/parse_multiple_chains.py --input_path={input_pdb} --output_path=$path_for_parsed_chains \

python3 $helperscripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "{binder_chain}" \

python3 $helperscripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "{binder_chain}" --position_list "$fixed_positions" \

python3 $runpath \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --path_to_model_weights $paramdir \
        --out_folder {out} \
        --num_seq_per_target {num_of_seq} \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
rm -r $TEMP
        """.format(
            input_pdb=rs_model02,
            out=sol_mpnn_output_dir,
            binder_chain=binder_chain,
            num_of_seq=seqs_to_generate
        )
        
        # make the slurm file
        with open(f"{slurm_output_folder}/{os.path.basename(rs_model02).rstrip('.pdb')}.slurm", 'w') as w01:
            w01.write(slurm_head)
            w01.write(slurm_hotspots)
            w01.write(slurm_body)
            w01.close()
            
    # submit slurm files
    for slurm_ in glob.iglob(os.path.join(slurm_output_folder, '*.slurm')):
        subprocess.call(['sbatch', slurm_])
       
    return None

## execution
if __name__ == '__main__':
    main()
