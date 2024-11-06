#!/bin/bash -l

#SBATCH --job-name=ipython-izar
#SBATCH --nodes 1
#SBATCH --gres gpu:1
#SBATCH --exclusive
#SBATCH --mem 20G
#SBATCH --time 10:00:00
#SBATCH --output jupyter-log-%J.out

module load gcc cuda cudnn mvapich2
source /home/khakzad/venv_sb/bin/activate

ipnport=$(shuf -i8000-9999 -n1)
 
jupyter-notebook --no-browser --port=${ipnport} --ip=$(hostname -i)
