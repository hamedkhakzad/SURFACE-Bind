#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 64G
#SBATCH --time 14:00:00
#SBATCH --account lpdi

module purge
module load gcc cmake python

source /work/lpdi/bin/load_masif_environment_jed.sh

PDBDIR="/work/lpdi/users/khakzad/Surfacome/pdbs/remained"


for FILE in $PDBDIR/*.pdb; do 
    n=${FILE##*/};
    m=${n%.pdb}
    NAME=$m"_A"
    
    ./data_prepare_one.sh --file $FILE $NAME;
    ./predict_site.sh $NAME;
    ./color_site.sh $NAME;
    ./compute_descriptors.sh $NAME;
    cp -r targets/template/ targets/$NAME
    cp -r targets/template_beta/ targets/$NAME"_beta"

done