BAM=$1
OUT=$2
CAPTURE=$3
UMI1=$4
UMI2=$5
IS=$6
LIBTYPE=$7
THREADS=$8
PREFIX=${BAM##*/}

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

if [[ $LIBTYPE == "DNA" ]]; then
  SCRIPT="filter_dna"
elif [[ $LIBTYPE == "RNA" ]]; then
  SCRIPT="filter_rna"
else
  echo "False library type argument"
  exit 1
fi

mkdir -p .tmp

samtools view -H $BAM > .tmp/"$PREFIX"_filtered
while IFS="" read -r p || [ -n "$p" ]
do
  FWD=$(cut -f2 <<< "$p")
  REV="`echo $(cut -f3 <<< "$p") | tr ACGTacgt TGCAtgca | rev`" 
  samtools view $BAM $(cut -f4 <<< "$p"):$(cut -f5 <<< "$p")-$(cut -f6 <<< "$p") | LC_ALL=C mawk -v fwd="$FWD" -v rev="$REV" -v umi1="$UMI1" -v umi2="$UMI2" -v is="$IS" -f $SCRIPTPATH/$SCRIPT.awk >> .tmp/"$PREFIX"_filtered
done < $CAPTURE

samtools view -@{threads} -Sb .tmp/"$PREFIX"_filtered > $OUT
rm .tmp/"$PREFIX"_filtered