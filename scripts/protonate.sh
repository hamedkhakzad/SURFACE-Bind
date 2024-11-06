ROOT="/work/lpdi/users/khakzad/Surfacome"
RAW="${ROOT}/pdbs/surfacome_error"
PROTONATED="${ROOT}/pdbs/surfacome_error_protonated"
REDUCE="/home/khakzad/bin/reduce"
LOG="${ROOT}/reduce.log"

mkdir -p $PROTONATED

IFS=$'\n'
set -f

COUNTER=0
TOTAL=$(ls $RAW | wc -l)

for PDB in $(ls $RAW);
do
  PROTEIN_PATH="${RAW}/${PDB}"
  PROTONATED_PROTEIN_PATH="${PROTONATED}/${PDB}"

  $REDUCE $PROTEIN_PATH > $PROTONATED_PROTEIN_PATH 2> $LOG

  COUNTER=$((COUNTER+1))
  PROGRESS=$((100 * COUNTER / TOTAL))
  BAR=$(seq -s= $PROGRESS|tr -d '[:digit:]')
  echo -ne "${BAR} ${COUNTER} (${PROGRESS}%)\r"
done

