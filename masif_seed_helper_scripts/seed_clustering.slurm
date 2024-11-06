#!/bin/bash -l
#SBATCH --job-name=seed_clustering
#SBATCH --nodes=1
#SBATCH --time=01:00:00
 
module load gcc/11.3.0 python/3.10.4
 
source /work/lpdi/users/khakzad/envs/jed_nb_3.6/bin/activate

site=$1

python3 -W ignore seed_clustering.py -s $site
 

