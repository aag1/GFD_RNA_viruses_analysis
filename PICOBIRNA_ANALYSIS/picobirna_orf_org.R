### R libraries
.libPaths('/groups/umcg-lld/tmp03/umcg-agulyaeva/R_LIB')
library('optparse')
sessionInfo()



### R functions
source('../function_plot_orf_organization.R')
plot_orf_organization



### input parameters
option_list = list(
    make_option('--protsF'),
	make_option('--rdrpF'),
    make_option('--orderF'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read data
rdrp <- read.table(
            opt$rdrpF,
            sep = '\t',
            header = TRUE,
            stringsAsFactors = FALSE)

order <- read.table(
            opt$orderF,
            header = FALSE,
            stringsAsFactors = FALSE)[, 1]

df <- read.table(
            opt$protsF,
            sep = '\t',
            header = FALSE,
            stringsAsFactors = FALSE)

colnames(df) <- c('orf_id', 'orf_from', 'orf_to', 'orf_strand', 'orf_info')

df$contig_id <- sub('_[0-9]+$', '', df$orf_id)

df$contig_len <- sub('^.+_length_([0-9]+)_.+$', '\\1', df$contig_id)
df$contig_len <- as.numeric(df$contig_len)

df$orf_func <- ifelse(df$orf_id %in% rdrp$Sequence, 'RdRp', '')



### plot ORF organization
pdf(
    'picobirna_orf_org.pdf',
    width = 16 / 2.54,
    height = 5 / 2.54
)

layout(matrix(1:15, ncol = 3))
par(cex = 0.5)


# contigs
for (x in order) {

    idx <- which(df$contig_id == x)
    if (any(df$orf_strand[idx] != df$orf_strand[idx][1])) { stop(paste(x, 'has ORFs on both strands!')) }


    tab <- data.frame(
                from = c(1, df$orf_from[idx]),
                to = c(df$contig_len[idx][1], df$orf_to[idx]),
                name = c('contig', df$orf_func[idx]),
                stringsAsFactors = FALSE
    )


    plot_orf_organization(tab, xmax = 2100)

    text(
        x = 0,
        y = 3.5,
        labels = paste(sub('_length_.+$', '', df$contig_id[idx][1]), ifelse(df$orf_strand[idx][1] == 1, '', 'rev. compl.')),
        adj = 0
    )

}


# axis
plot(
    NA,
    xlim = c(-0.1, 1.1) * 2100,
    ylim = c(-0.5, 4.5),
    xaxs = 'i', yaxs = 'i', 
    axes = FALSE, ann = FALSE, bty = 'n'
)

axis(
    side = 1,
    pos = 3.5,
    at = seq(0, 2000, 400),
    mgp = c(3, 0.5, 0)
)

text(
    x = 1000,
    y = 0.5,
    labels = 'Contig length, nt'
)


dev.off()
