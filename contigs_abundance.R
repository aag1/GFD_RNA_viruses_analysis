### R libraries
.libPaths('/groups/umcg-lld/tmp03/umcg-agulyaeva/R_LIB')
library('optparse')
sessionInfo()



### R functions
source('function_plot_contigs_abundance.R')
plot_contigs_abundance



### input parameters
option_list = list(
	make_option('--rdrp_contigs'),
	make_option('--aichi_contigs'),
	make_option('--pbv_contigs'),
	make_option('--counts_table'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read input
r <- read.table(
				opt$rdrp_contigs,
				sep = '\t',
				header = TRUE,
				stringsAsFactors = FALSE
)

a <- read.table(
				opt$aichi_contigs,
				sep = '\t',
				header = TRUE,
				stringsAsFactors = FALSE
)

p <- read.table(
				opt$pbv_contigs,
				header = FALSE,
				stringsAsFactors = FALSE
)[,1]

tab <- read.table(
				opt$counts_table,
				sep = '\t',
				row.names = 1,
				header = TRUE,
				stringsAsFactors = FALSE
)



### function
get_abundance <- function (sele, tab) {

	d <- tab[sele,]
	d <- t(as.matrix(d))
	return(d)

}



### calculate contigs abundance
L <- list(
	aichi_virus = unique(a$qseqid),
	picobirna = p,
	tombus_like = r$contig_id[grep('tombus-like', r$RefSeq_prot_taxo)],
	plant = r$contig_id[grep('Virgaviridae|Tymovirales|Luteoviridae', r$RefSeq_prot_taxo)]
)


for (x in names(L)) {
	
	print(all(L[[x]] %in% rownames(tab)))

	d <- get_abundance(L[[x]], tab)

	write.table(
		d,
		sep = '\t',
		quote = FALSE,
		file = paste0(x, '_contigs_abundance.txt')
	)

}



### plot contigs abundance
pdf(
    paste0('picobirna_tombus_contigs_abundance.pdf'),
    width = 18 / 2.54,
    height = 20 / 2.54
)

layout(matrix(1:2, nrow = 1), width = c(14, 4))

par(ps = 10)


for (x in c('picobirna', 'tombus_like')) {

    if (x=='picobirna')   {
                            par(mar = c(2, 3, 8, 3))
                            BR <- c(0, 10^seq(1, 3, 0.06))
                            BRL <- c(0, 10^seq(1, 3, 1))
    }
    if (x=='tombus_like') {
                            par(mar = c(2, 0, 8, 1))
                            BR <- c(0, 10^seq(1, 5, 0.1))
                            BRL <- c(0, 10^seq(1, 5, 1))
    }


	d <- get_abundance(L[[x]], tab)

    rownames(d) <- sub('^GFD_', '', rownames(d))
    d <- d[order(as.numeric(rownames(d))),]

    colnames(d) <- sub('^(.+)_length_.+$', '\\1', colnames(d))


    plot_contigs_abundance(
			d,
			sample_space = 3,
			brackets = BR,
            brackets_lab = BRL
    )

}


par(new=TRUE, mar=rep(0, 4))
layout(matrix(1))
plot(NA, xlim=c(0, 18), ylim=c(0, 1), xaxs='i', yaxs='i', axes=FALSE, ann=FALSE)
text(x=c(0.75, 13.5), y=0.95, labels=c('A', 'B'), cex=2.25)


dev.off()
