# GFD RNA viruses analysis
This repository contains scripts used to analyse RNA viruses of the GFD virome project:

* __run_hmmer.sh__ uses HMMER to compare proteins of all non-redundant contigs to PFAM RdRp profiles

* __process_hmmer_out.sh__ and __process_hmmer_out.R__ collect hits characterized by E-value <0.001 and with RdRp motifs A-C into a single table, place corresponding alignments into a single file

* __run_blastp.sh__ uses BLASTP to compare identified RdRp proteins to Viral RefSeq proteins

* __process_blastp_out.sh__ retrieves taxonomy of the viruses that provided hits on the previous step by __get_RefSeq_proteins_info.pl__, and summarizes results obtained on the previous step by __process_blastp_out.R__

* __contigs_abundance.sh__ launches __contigs_abundance.R__ that (1) retrieves RPKM counts for the four groups of RdRp-containing contigs and (2) plots RPKM counts for picobirnavirus and tombus-like contigs

* folder __PICOBIRNA_ANALYSIS/__:
    * __build_picobirna_msa_tree.sh__ builds picobirnavirus RdRp MSA and tree
    * __nexus2newick.sh__ converts tree file from NEXUS to NEWICK format following midpoint-rooting by FigTree
    * __plot_picobirna_tree.R__ plots tree and prepares an MSA of representatives for plotting
    * __picobirna_orf_org.sh__ launches __picobirna_orf_org.R__ that plots ORF organization of the picobirnavirus contigs

* folder __TOMBUSLIKE_ANALYSIS/__:
    * __build_tombus_msa.sh__ builds tombus-like RdRp MSA
    * __tombus_orf_org.sh__ launches __tombus_orf_org.R__ that plots ORF organization of the tombus-like contigs
