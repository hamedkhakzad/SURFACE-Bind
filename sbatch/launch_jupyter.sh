#!/bin/bash
#SBATCH --job-name=ipython-notebook
#SBATCH --nodes=1
#SBATCH --cpus-per-task 4
#SBATCH --mem 32G
#SBATCH --time=09:00:00
#SBATCH --output jupyter-log-%J.out
 
module load gcc/8.4.0 mvapich2
module load py-torch/1.6.0-openmp
 
source ~/venv_fidis_nb/bin/activate
 
ipnport=$(shuf -i8000-9999 -n1)
 
jupyter-notebook --no-browser --port=${ipnport} --ip=$(hostname -i)
