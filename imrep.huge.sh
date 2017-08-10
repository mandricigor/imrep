
##########################################################################################
####    This is the wrapper script for ImReP to handle huge fasta/fastq read files    ####
##########################################################################################


inputfile=$1
nrchunks=$2
outputfile=$3

DATE=`date +%Y-%m-%d-%H-%M-%S`

tmpdir=/tmp/experiment-${DATE}
inputs=$tmpdir/inputs
outputs=$tmpdir/outputs

mkdir -p $tmpdir
mkdir -p $inputs
mkdir -p $outputs

prefix=$inputs/tmpfile
split -n $nrchunks -d -e $inputfile $prefix


for f in `ls $inputs`; do 
    echo "python imrep.py --noCast --fastq $inputs/$f $outputs/$f &" >> $tmpdir/runAll.sh
done
echo "wait" >> $tmpdir/runAll.sh


bash $tmpdir/runAll.sh > /dev/null


for f in `ls $outputs`; do
    tail -n +2 $outputs/$f >> $tmpdir/merged.out
done


python cast.py $tmpdir/merged.out $outputfile

rm -rf $tmpdir

echo "DONE"


