#!/bin/bash -l
#SBATCH --job-name=ipython-notebook
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --output jupyter-log-%J.out

module purge
module load gcc python
 
source /work/lpdi/users/khakzad/envs/jed_nb_2/bin/activate
 
ipnport=$(shuf -i8000-9999 -n1)
 
jupyter-notebook --no-browser --port=${ipnport} --ip=$(hostname -i)
