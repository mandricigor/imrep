#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1

parser.add_argument('bam',help='bam with mapped and unmapped reads, BAM file needs to be indexed')

parser.add_argument('-s', '--species', default='human', type=str,
help='species (human or mouse, default human)')
parser.add_argument('-r', '--release', default='hg19', type=str,
help='For human we support hg19 (default) and hg38. For mouse we support only mm10 (default)')
parser.add_argument('-chr', '--chr', default='', type=str,
help='Use this option to specify the format of chromosome name in the bam file. By default it is : 1,2,3,...,X. In case the format is : chr1, chr2,..,chrX, use : -chr chr)')
parser.add_argument('-c', '--chains', default='IGH,IGK,IGL,TRA,TRB,TRD,TRG', type=str,
help='chains: comma separated values from IGH,IGK,IGL,TRA,TRB,TRD,TRG')

EOF




#echo ${SPECIES[@]}
#echo required infile: $BAM
#echo required basenamefile: $basenamePUT_CLONES



bam=$BAM


basename=$(echo  ${bam##*/} | awk -F ".bam" '{print $1}')
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
chr=${CHR[@]}

mkdir $basename

samtools=${DIR}/ext/samtools


if [ ! -f ${bam}.bai ]; then
echo "Warning : Index of $bam not found"
echo "Creating index ... (samtools index $bam)"

samtools index $bam
echo "DONE!"

fi



#--------------------------------------------------------------------------------------
#UNMAPPED READS
$samtools view -f 0x4 -bh  $bam | $samtools bam2fq - > ${basename}/${basename}_unmapped.fastq



#--------------------------------------------------------------------------------------
#MAPPED READS

#human - hg38
#extract reads from IGH chr14, 105586437..106879844
#extract reads from IGK chr2, 88857361..90235368
#extract reads from IGL chr22, 22026076..22922913
#extract reads from TCRA chr14, 21621904..22552132
#extract reads from TCRB chr7, 142299011..142813287basename
# reads from TCRD inside TCRA
# #extract reads from TCRG chr7, 38240024..38368055


##mm10
## TCRG : 13	NC_000079.6 (19178042..19356476)

if [ ${SPECIES[@]} == "human" ]
then
echo "-->Human BCR and TCR annotations are used"


if [ ${RELEASE[@]} == "hg19" ]
then
echo "-->Release hg19"




if [[ ${CHAINS[@]} == *"IGH"* ]]
then
$samtools view -bh ${bam} ${chr}14:106032614-107288051 | $samtools view -bh -F 4 - | $samtools bam2fq - >${basename}/${basename}_mapped_immune.fastq
fi

if [[ ${CHAINS[@]} == *"IGK"* ]]
then
$samtools view -bh ${bam} ${chr}2:89156874-89630436 | $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
fi

if [[ ${CHAINS[@]} == *"IGL"* ]]
then
$samtools view -bh ${bam} ${chr}22:22380474-23265085 | $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
fi

if [[ ${CHAINS[@]} == *"TRA"* ]]
then
$samtools view -bh ${bam} ${chr}14:22090057-23021075 | $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq

fi

if [[ ${CHAINS[@]} == *"TRB"* ]]
then
$samtools view -bh ${bam} ${chr}7:141998851-142510972 | $samtools view -bh -F 4 - | $samtools bam2fq -    >>${basename}/${basename}_mapped_immune.fastq
fi

if [[ ${CHAINS[@]} == *"TRG"* ]]
then
$samtools view -bh ${bam} ${chr}7:38279625-38407656 | $samtools view -bh -F 4 - | $samtools bam2fq - >>${basename}/${basename}_mapped_immune.fastq
fi

fi

if [ ${RELEASE[@]} == "hg38" ]
then
echo "-->Release hg38"

$samtools view -bh ${bam} ${chr}14:105586437-106879844 | $samtools view -bh -F 4 - | $samtools bam2fq - >${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}2:88857361-90235368 | $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}22:22026076-22922913 | $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}14:21621904-22552132 | $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}7:142299011-1428132872 | $samtools view -bh -F 4 - | $samtools bam2fq -    >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}7:38240024-38368055 | $samtools view -bh -F 4 - | $samtools bam2fq - >>${basename}/${basename}_mapped_immune.fastq
fi


else
echo "-->Mouse(mm10) BCR and TCR annotations are used"
#igh="chr12:113258768-116009954"
#igk="chr6:67555636-70726754"
#igl="chr16:19026858-19260844"
#tcra="chr14:52427967-54224198"
#tcrb="chr6:40891296-41558371"
#tcrg="chr13:19178042-19356476"
#tcrd is inside tcra

$samtools view -bh ${bam} ${chr}12:113258768-116009954 | $samtools view -bh -F 4 - | $samtools bam2fq - >${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}6:67555636-70726754| $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}16:19026858-19260844| $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}14:52427967-54224198 | $samtools view -bh -F 4 - | $samtools bam2fq -  >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}6:40891296-41558371 | $samtools view -bh -F 4 - | $samtools bam2fq -    >>${basename}/${basename}_mapped_immune.fastq
$samtools view -bh ${bam} ${chr}13:19178042-19356476 | $samtools view -bh -F 4 - | $samtools bam2fq - >>${basename}/${basename}_mapped_immune.fastq
fi



n=$(wc -l ${basename}/${basename}_mapped_immune.fastq | awk '{print $1/4}')
echo "Number of mapped reads extracted from BCR and TCR loci : "$n
n_unm=$(wc -l ${basename}/${basename}_unmapped.fastq | awk '{print $1/4}')
echo "Number of unmapped reads extracted: "$n_unm

N=$((n + n_unm))

echo "Merging mapped reads from BCR/TCR loci and unmapped reads ..."
echo "In total $N reads were available for ImReP analysis"




cat ${basename}/${basename}_mapped_immune.fastq ${basename}/${basename}_unmapped.fastq > ${basename}/${basename}_unmapped_plus_immune.fastq


rm ${basename}/${basename}_mapped_immune.fastq ${basename}/${basename}_unmapped.fastq


if [ ${CHAINS[@]} == "IGH,IGK,IGL,TRA,TRB,TRD,TRG" ]
then
echo "Profile all chains : IGH,IGK,IGL,TRA,TRB,TRD,TRG"
python ${DIR}/imrep.py --fastq ${basename}/${basename}_unmapped_plus_immune.fastq ${basename}/${basename}.cdr3
else
echo "Profile only chains :" ${CHAINS[@]}

python ${DIR}/imrep.py -c ${CHAINS[@]} --fastq ${basename}/${basename}_unmapped_plus_immune.fastq ${basename}/${basename}.cdr3
fi



echo "Results are in ${basename}"
echo "DONE!"

