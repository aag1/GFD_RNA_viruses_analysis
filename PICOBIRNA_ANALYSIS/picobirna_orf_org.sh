#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2.out
#SBATCH --error=job2.err
#SBATCH --mem=2gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### load modules
module load seqtk
module load R/3.5.1-foss-2015b-bare
module list



### get ORF organization info
allF=../nonredundant_contigs.min1000.AA.fasta

idsF=GFD_picobirna_RdRp_contigs_ids.txt
grep 'Picobirnaviridae' ../GFD_contigs_with_RdRp.txt |
cut -d$'\t' -f2 > ${idsF}

protsF=GFD_picobirna_RdRp_contigs_proteins.txt
grep -F -f ${idsF} ${allF} |
sed 's/^>//' |
sed 's/ \# /\t/g' > ${protsF}



### plot ORF organization
Rscript picobirna_orf_org.R \
    --protsF ${protsF} \
    --rdrpF ../GFD_contigs_with_RdRp.txt \
    --orderF picobirna_contigs_tree_order.txt
