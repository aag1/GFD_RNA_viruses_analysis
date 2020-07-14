#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1.out
#SBATCH --error=job1.err
#SBATCH --mem=5gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=10
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### load modules
module load seqtk
module load Anaconda3
module list



### get GFD picobirna RdRps
allF=../nonredundant_contigs.min1000.AA.fasta

idsF=GFD_picobirna_RdRp_ids.txt
grep 'Picobirnaviridae' ../GFD_contigs_with_RdRp.txt |
cut -d$'\t' -f1 |
sed 's/$/ /' > ${idsF}

long_idsF=GFD_picobirna_RdRp_ids_long.txt
grep -F -f ${idsF} ${allF} |
sed 's/^>//' > ${long_idsF}

seleF='GFD_picobirna_RdRp_seq.fasta'
seqtk subseq -l 60 ${allF} ${long_idsF} > ${seleF}



### get outgroup sequence
outgrF='picobirna_RdRp_outgroup.fasta'
if [ -f ${outgrF} ]; then rm ${outgrF}; fi

wget -O ${outgrF} 'https://www.ncbi.nlm.nih.gov/search/api/sequence/APG78210.1/?report=fasta'
sed -i '/^\s*$/d' ${outgrF}
sed -i 's/^>.\+$/>Hubei_earwig_virus_2/' ${outgrF}



### add GFD picobirna RdRps to the ICTV ones
source activate '/groups/umcg-lld/tmp03/umcg-agulyaeva/CONDA/envs/MAFFT'
conda list

ictv_msa='/groups/umcg-lld/tmp03/umcg-agulyaeva/DATA/Picobirnaviridae_RdRp_alignment_ICTV2018/ODR.Picobirna.Fig4.RdRP.align_NO_PARTITIVIRUS.fas'
join_msa='ICTV_GFD_picobirna_RdRp_msa.fasta'
mafft --thread 10 --add ${seleF} ${ictv_msa} > ${join_msa}



### add outgroup
join_msa2='ICTV_GFD_picobirna_RdRp_msa_outgr.fasta'
mafft --thread 10 --add ${outgrF} ${join_msa} > ${join_msa2}

conda deactivate



### build tree
source activate '/groups/umcg-lld/tmp03/umcg-agulyaeva/CONDA/envs/IQTREE'
conda list

tree_file='ICTV_GFD_picobirna_RdRp_tree.nwk'
iqtree -s ${join_msa2} -bb 1000 -nt 10

conda deactivate
