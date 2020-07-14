#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1.out
#SBATCH --error=job1.err
#SBATCH --mem=5gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=5
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### load modules
module load seqtk
module load Anaconda3
module list



### file for unaligned sequences
seqF='tombus_like_RdRp_seq.fasta'
if [ -f ${seqF} ]; then rm ${seqF}; fi



### get RefSeq tombus-like proteins
sele=( $(grep 'tombus-like' ../GFD_contigs_with_RdRp.txt | cut -d$'\t' -f8) )

for id in ${sele[@]}
do
	wget -O ${id}'.fasta' 'https://www.ncbi.nlm.nih.gov/search/api/sequence/'${id}'/?report=fasta'
	cat ${id}'.fasta' | sed '/^\s*$/d' >> ${seqF}
done



### Get GFD tombus-like proteins
allF='../nonredundant_contigs.min1000.AA.fasta'

idsF='GFD_tombus_like_RdRp_ids.txt'
grep 'tombus-like' ../GFD_contigs_with_RdRp.txt |
cut -d$'\t' -f1 |
sed 's/$/ /' > ${idsF}

long_idsF='GFD_tombus_like_RdRp_ids_long.txt'
grep -F -f ${idsF} ${allF} |
sed 's/^>//' > ${long_idsF}

seqtk subseq -l 60 ${allF} ${long_idsF} >> ${seqF}



### build MSA
source activate '/groups/umcg-lld/tmp03/umcg-agulyaeva/CONDA/envs/MAFFT'
conda list

mafft --thread 5 ${seqF} > tombus_like_RdRp_msa.fasta

conda deactivate



### rename sequences in MSA
perl -i -pe 's/^>(YP_009342273.1).*$/>WTLV17_\1/' tombus_like_RdRp_msa.fasta
perl -i -pe 's/^>(YP_009337712.1).*$/>HTLV36_\1/' tombus_like_RdRp_msa.fasta
perl -i -pe 's/^(>GFD_.+)_length_.+$/\1/' tombus_like_RdRp_msa.fasta
