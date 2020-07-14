#!/bin/bash
#SBATCH --job-name=hmmer
#SBATCH --output=hmmer.out
#SBATCH --error=hmmer.err
#SBATCH --mem=10gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=10
#SBATCH --export=NONE
#SBATCH --get-user-env=L


# file with protein content of non-redundant contigs
cp -p /groups/umcg-lld/tmp03/umcg-sgarmaeva/GFD/GFD_all/Roux_analysis/pVOGs/nonredundant_contigs.min1000.AA.fasta .
db_aa='nonredundant_contigs.min1000.AA.fasta'


# create folder for HMMER output
dir='hmmer_out/'
if [ -d ${dir} ]; then rm -rf ${dir}; fi
mkdir ${dir}
cd ${dir}


# run HMMER
module load Anaconda3
module list
source activate '/groups/umcg-lld/tmp03/umcg-agulyaeva/CONDA/envs/HMMER'
conda list


for f in ../RdRp_MSAs/*fasta
do
	b=$(basename $f '.fasta')
	h=${b}'.hmm'
	o=${b}'_vs_GFD.txt'
	t=${b}'_vs_GFD_table.txt'

	hmmbuild ${h} ${f}
	hmmsearch --cpu 10 ${h} ../${db_aa} > ${o}

	sed -n '/--- full sequence ---/,/------ inclusion threshold ------/p' ${o} |
	sed '1d; 3d; $d' |
	sed 's/^ \+//' |
	sed 's/  \+/\t/g' > ${t}
done


callanan2020_hmm='/groups/umcg-lld/tmp03/umcg-agulyaeva/DATA/Callanan2020_PMID32083183_SupplData/Supplementary_data_1/ssRNAphage_detection_HMM5-MC/hmm_m5-mc'
hmmsearch --cpu 10 ${callanan2020_hmm} ../${db_aa} > 'callanan2020_hmm_vs_GFD.txt'


conda deactivate
