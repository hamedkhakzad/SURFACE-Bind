#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 20G
#SBATCH --time 10:00:00
#SBATCH --account lpdi

module purge
module load gcc cmake python
source /home/khakzad/venv_sb/bin/activate
cd /work/lpdi/users/khakzad/Surfacome

bash scripts/protonate.sh

