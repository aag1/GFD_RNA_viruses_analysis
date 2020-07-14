#!/bin/bash
#SBATCH --job-name=blastp
#SBATCH --output=blastp.out
#SBATCH --error=blastp.err
#SBATCH --mem=5gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=10
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### load software packages
module load seqtk
module load BLAST+
module list



### create folder for blastp output
dir=blastp_out/
if [ -d ${dir} ]; then rm -rf ${dir}; fi
mkdir ${dir}
cd ${dir}



### get proteins of interest
allF=../nonredundant_contigs.min1000.AA.fasta

idsF=GFD_proteins_with_HMMER_RdRp_hits_ids.txt
cut -d$'\t' -f1 ../GFD_contigs_with_HMMER_RdRp_hits.txt |
sed '1d' |
sed 's/$/ /' > ${idsF}

long_idsF=GFD_proteins_with_HMMER_RdRp_hits_ids_long.txt
grep -F -f ${idsF} ${allF} |
sed 's/^>//' > ${long_idsF}

seleF=GFD_proteins_with_HMMER_RdRp_hits.fasta
seqtk subseq -l 60 ${allF} ${long_idsF} > ${seleF}



### run tblastx
outF=GFD_proteins_with_HMMER_RdRp_hits_VS_VirProtRefSeq98.txt
blastp \
	-query ${seleF} \
	-db /groups/umcg-lld/tmp03/umcg-agulyaeva/DATA/refseq98_vir_prot/refseq98_vir_prot.fasta \
	-evalue 1e-3 \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
	-num_threads 10 \
	-out ${outF}
