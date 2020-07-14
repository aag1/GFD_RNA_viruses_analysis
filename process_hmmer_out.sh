#!/bin/bash
#SBATCH --job-name=jobA
#SBATCH --output=jobA.out
#SBATCH --error=jobA.err
#SBATCH --mem=5gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### load modules
module load R/3.5.1-foss-2015b-bare
module list



# select HMMER hits
Rscript process_hmmer_out.R \
	--hmmer_out_dir 'hmmer_out/' \
	--vir_cont_orig 'from_Resilio/Contigs_metadata.txt'



# extract corresponding HMMER alignments
fileO='GFD_contigs_with_HMMER_RdRp_hits.ali'
echo -e '\n\n\n' > ${fileO}
sed '1d' GFD_contigs_with_HMMER_RdRp_hits.txt | while read -r line
do
	regexp=$(echo ${line} | cut -d' ' -f1)' '

	rdrp=$(echo ${line} | cut -d' ' -f6)

	fileI=$(ls hmmer_out/Pfam32.0_${rdrp}_*_full_vs_GFD.txt)

	sed -n "/^>> ${regexp}/,/^>> /p" ${fileI} | sed '$d' >> ${fileO}

	echo -e '\n\n\n' >> ${fileO}

done



# remove hit w/o catalityc motif C
sed -i '/^GFD_14\.3_NODE_17690_length_1579_cov_4\.381234_1\t/d' GFD_contigs_with_HMMER_RdRp_hits.txt
