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
	make_option('--rdrpF'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read data
rdrp <- read.table(
            opt$rdrpF,
            sep = '\t',
            header = TRUE,
            stringsAsFactors = FALSE)

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
    'tombus_like_orf_org.pdf',
    width = 16 / 2.54,
    height = 5 / 2.54
)

layout(matrix(1:3, ncol = 1))


# contigs
for (x in unique(df$contig_id)) {

    idx <- which(df$contig_id == x)
    if (any(df$orf_strand[idx] != df$orf_strand[idx][1])) { stop(paste(x, 'has ORFs on both strands!')) }


    tab <- data.frame(
                from = c(1, df$orf_from[idx]),
                to = c(df$contig_len[idx][1], df$orf_to[idx]),
                name = c('contig', df$orf_func[idx]),
                stringsAsFactors = FALSE
    )


    plot_orf_organization(tab, xmax = 3000)

    if (x == 'GFD_11.7_NODE_3204_length_2878_cov_147.917109') {
        rect(
            xleft = 201,
            xright = 548,
            ybottom = 0,
            ytop = 1,
            col = 'grey90',
            density = 20,
            border = FALSE
        )
        rect(
            xleft = 201,
            xright = 548,
            ybottom = 0,
            ytop = 1
        )
    }
    if (x == 'GFD_14.1_NODE_11134_length_2601_cov_22.109976') {
        rect(
            xleft = 1946,
            xright = 2176,
            ybottom = 0,
            ytop = 1,
            col = 'grey90',
            density = 20,
            border = FALSE
        )
        rect(
            xleft = 1946,
            xright = 2176,
            ybottom = 0,
            ytop = 1
        )
    }

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
    xlim = c(-0.1, 1.1) * 3000,
    ylim = c(-0.5, 4.5),
    xaxs = 'i', yaxs = 'i', 
    axes = FALSE, ann = FALSE, bty = 'n'
)

axis(
    side = 1,
    pos = 3.5,
    at = seq(0, 3000, 500),
    mgp = c(3, 1, 0)
)

text(
    x = 1500,
    y = 0.5,
    labels = 'Contig length, nt'
)


dev.off()
