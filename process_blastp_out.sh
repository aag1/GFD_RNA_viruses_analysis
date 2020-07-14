#!/bin/bash
#SBATCH --job-name=jobB
#SBATCH --output=jobB.out
#SBATCH --error=jobB.err
#SBATCH --mem=5gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### load modules (order of Perl modules matter!)
module load Perl
module load BioPerl
module load PerlPlus
module load R/3.5.1-foss-2015b-bare
module list



### array of unique RefSeq protein ids
blastp_out_file='blastp_out/GFD_proteins_with_HMMER_RdRp_hits_VS_VirProtRefSeq98.txt'
refseq_prot_ids=( $(cut -d$'\t' -f2 ${blastp_out_file} | sort | uniq) )



### retrieve info about RefSeq proteins
perl get_RefSeq_proteins_info.pl ${refseq_prot_ids[@]} > 'VirProtRefSeq98_info.txt'



### prepare final table
Rscript process_blastp_out.R \
	--hmmer_table 'GFD_contigs_with_HMMER_RdRp_hits.txt' \
	--blastp_out ${blastp_out_file} \
	--refseq_info 'VirProtRefSeq98_info.txt'
