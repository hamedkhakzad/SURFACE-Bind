#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --gres gpu:1
#SBATCH --mem 64G
#SBATCH --time 60:00:00
#SBATCH --account lpdi

module purge
module load gcc cmake python cuda cudnn
source /home/khakzad/venv_sb/bin/activate
cd /work/lpdi/users/khakzad/Surfacome

python -W ignore prepare_dataset.py \
                 --data /work/lpdi/users/khakzad/Surfacome/database \
                 --structures /work/lpdi/users/khakzad/Surfacome/pdbs/surfacome_chopped_protonated \
                 --dmasif_model /work/lpdi/users/khakzad/dMaSIF/models/dMaSIF_site_3layer_16dims_9A_100sup_epoch64 \
                 --experiment_name None \
                 --embedding_layer dMaSIF \
                 --site True \
                 --single_protein True \
                 --n_layers 3 \
                 --emb_dims 16 \
                 --sup_sampling 100 \
                 --resolution 1 \
                 --device cuda:0

