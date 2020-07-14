### R libraries
.libPaths('/groups/umcg-lld/tmp03/umcg-agulyaeva/R_LIB')
library('optparse')
sessionInfo()




### input parameters
option_list = list(
	make_option('--hmmer_table'),
	make_option('--blastp_out'),
	make_option('--refseq_info'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)




### read input
DF <- read.table(
				opt$hmmer_table,
				sep = '\t',
				header = TRUE,
				stringsAsFactors = FALSE
)


bp <- read.table(
				opt$blastp_out,
				sep = '\t',
				header = FALSE,
				stringsAsFactors = FALSE
)
colnames(bp) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen')
bp$qcov <- (bp$qend - bp$qstart + 1) / bp$qlen * 100


rs <- read.table(
				opt$refseq_info,
				sep = '\t',
				header = TRUE,
				stringsAsFactors = FALSE
)




### add columns to the table
DF$RefSeq_prot_id <- NA
DF$RefSeq_prot_eval <- NA
DF$RefSeq_prot_pident <- NA
DF$RefSeq_prot_qcov <- NA
DF$RefSeq_prot_len <- NA
DF$RefSeq_prot_desc <- NA
DF$RefSeq_prot_taxo <- NA

colnam1 <- c('RefSeq_prot_id', 'RefSeq_prot_eval', 'RefSeq_prot_pident', 'RefSeq_prot_qcov')
colnam2 <- c('sseqid', 'evalue', 'pident', 'qcov')
colnam3 <- c('RefSeq_prot_len', 'RefSeq_prot_desc', 'RefSeq_prot_taxo')

for (i in 1:nrow(DF)) {

	# select BLASTP hit
	idx <- which(bp$qseqid == DF$Sequence[i] & bp$qcov > 75)
	if (length(idx) == 0) { next }

	sele <- idx[ bp$pident[idx] == max(bp$pident[idx]) ][1]


	# fill in info about BLASTP hit
	DF[i, colnam1] <- bp[sele, colnam2]


	# fill in info about RefSeq target
	n <- which(rs$RefSeq_prot_id == DF$RefSeq_prot_id[i])
	if (length(n) == 1) {

		DF[i, colnam3] <- rs[n, colnam3]
		DF$RefSeq_prot_desc <- sub(' \\[.+\\]\\.$', '', DF$RefSeq_prot_desc)

	}

}

DF$RefSeq_prot_pident <- round(DF$RefSeq_prot_pident, 1)
DF$RefSeq_prot_qcov <- round(DF$RefSeq_prot_qcov, 1)




### Write table
write.table(
	DF,
	sep = '\t',
	row.names = FALSE,
	quote = FALSE,
	file = 'GFD_contigs_with_RdRp.txt'
)
