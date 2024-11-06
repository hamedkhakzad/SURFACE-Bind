#!/bin/bash
#SBATCH --job-name=surfacome
#SBATCH --nodes=1
#SBATCH --cpus-per-task 4
#SBATCH --mem 32G
#SBATCH --time=36:00:00

module load gcc/11.3.0 python/3.10.4

source /work/lpdi/users/khakzad/envs/jed_nb_3.6/bin/activate

cd /work/lpdi/users/khakzad/Surfacome

python -W ignore analyze_dataset.py \
                 --data /work/lpdi/users/khakzad/Surfacome/database \
                 --structures /work/lpdi/users/khakzad/Surfacome/pdbs/target_list_1_chopped_protonated \
                 --outpath /work/lpdi/users/khakzad/Surfacome/outputs_117 \
                 --threshold 0.7 \
                 --noPointscluster 200 \
                 --device cpu

