### R libraries
.libPaths('/groups/umcg-lld/tmp03/umcg-agulyaeva/R_LIB')
library('optparse')
sessionInfo()




### input parameters
option_list = list(
	make_option('--hmmer_out_dir'),
	make_option('--vir_cont_orig'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)




### read input
orig <- read.table(
				opt$vir_cont_orig,
				sep = '\t',
				header = TRUE,
				row.names = 1,
				stringsAsFactors=FALSE
)
orig <- orig[, c('Refseq', 'pVOGs', 'Circ', 'VirSorter', 'crass', 'dark_matter')]


hmmer_out_files <- list.files(
						path = opt$hmmer_out_dir,
						pattern = '_table.txt$',
						full.names = TRUE
)

DF <- NULL

for (f in hmmer_out_files) {

	df <- read.table(
				f,
				sep = '\t',
				header = TRUE,
				comment.char = '',
				stringsAsFactors=FALSE
	)
	colnames(df)[colnames(df) == 'E.value'] <- 'rdrp_eval'
	if (nrow(df) == 0) {next}


	df$rdrp_name <- gsub('(^Pfam32.0_)|(_PF[0-9]{5}_full_vs_GFD_table.txt$)', '', basename(f))


	# filter out hits with E-value >= 0.001
	df <- df[df$rdrp_eval < 0.001, ]


	# parse df$Description column
	# https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output#gene-coordinates
	l <- strsplit(paste0(' ', df$Description), ' # ')
	df$orf_coo <- unlist(lapply(l, function (v) { paste0(v[2], '-', v[3], ';', ifelse(as.numeric(v[4]) == 1, 'forward', 'reverse')) }))


	# filter out ORFs <= 50 aa
	df$orf_aa <- unlist(lapply(l, function (v) { (diff(as.numeric(v[2:3])) + 1) / 3 }))
	df <- df[df$orf_aa > 50, ]


	DF <- rbind(DF, df)

}

DF$contig_id <- sub('_[0-9]+$', '', DF$Sequence)




### For each protein, preserve info about RdRp hit(s) with the lowest E-value
for (s in unique(DF$Sequence)) {

	idx <- which(DF$Sequence == s)

	if (length(idx) == 1) {next}

	E <- DF$rdrp_eval[idx]
	N <- which(E == min(E))
	bad <- idx[ -N ]
	DF <- DF[-bad, ]

}
cat(
	'\nA single hit was preserved for each protein in the table:',
	nrow(DF) == length(unique(DF$Sequence)),
	'\n\n'
)




### Add info about contig fate in GFD pipeline
DF$contig_fate <- 'Discarded'
idx <- which(DF$contig_id %in% rownames(orig))
DF$contig_fate[idx] <- sapply(DF$contig_id[idx], function (x) {
								src <- colnames(orig)[ orig[x,] != 0 ]
								src <- paste(src, collapse = ';')
								return(src)
})




### Write table
write.table(
	DF[, c('Sequence', 'contig_id', 'contig_fate', 'orf_coo', 'orf_aa', 'rdrp_name', 'rdrp_eval')],
	sep = '\t',
	row.names = FALSE,
	quote = FALSE,
	file = 'GFD_contigs_with_HMMER_RdRp_hits.txt'
)
