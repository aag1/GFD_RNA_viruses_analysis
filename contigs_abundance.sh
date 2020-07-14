#!/bin/bash
#SBATCH --job-name=jobC
#SBATCH --output=jobC.out
#SBATCH --error=jobC.err
#SBATCH --mem=1gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --qos=dev



### load modules
module load R/3.5.1-foss-2015b-bare
module list



### Aichi virus contigs
# NC_001918.1 3Dpol 6608 - 8011 nt
# NC_004421.1 3Dpol 6791 - 8197 nt
echo -e 'qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tqstart\tqend\tsstart\tsend\tstitle' > Aichi_virus_contigs.txt
grep 'Aichi' RNA_VIRUS_CONTIGS_IDENTIFIED_BY_REFSEQ_BUT_NOT_RDRP/nr_contigs_viral_refseq_outfmt6.0819.txt >> Aichi_virus_contigs.txt



### RdRp contigs abundance
Rscript contigs_abundance.R \
	--rdrp_contigs 'GFD_contigs_with_RdRp.txt' \
	--aichi_contigs 'Aichi_virus_contigs.txt' \
	--pbv_contigs 'PICOBIRNA_ANALYSIS/picobirna_contigs_tree_order.txt' \
	--counts_table 'from_Resilio/RPKM_counts_GFD.txt'
