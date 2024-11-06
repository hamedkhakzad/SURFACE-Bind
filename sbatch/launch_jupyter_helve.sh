#!/bin/bash
#SBATCH --job-name=ipython-notebook
#SBATCH --nodes=1
#SBATCH --cpus-per-task 4
#SBATCH --mem 32G
#SBATCH --time=5:00:00
#SBATCH --output jupyter-log-%J.out
 
module load gcc mvapich2
 
source /work/lpdi/users/khakzad/envs/helve_nb/bin/activate
 
ipnport=$(shuf -i8000-9999 -n1)
 
jupyter-notebook --no-browser --port=${ipnport} --ip=$(hostname -i)
