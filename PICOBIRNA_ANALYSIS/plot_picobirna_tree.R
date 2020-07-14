TREE_FILE <- 'ICTV_GFD_picobirna_RdRp_msa_outgr.FigTreeMidpoint.nwk'
MSA_FILE <- 'ICTV_GFD_picobirna_RdRp_msa_outgr.fasta'



### R libraries
.libPaths('/groups/umcg-lld/tmp03/umcg-agulyaeva/R_LIB')
library('ape')
library('seqinr')
sessionInfo()



### R functions
source('../function_tree_support_pch.R')
tree_support_pch
source('../function_remove_empty_columns.R')
remove_empty_columns



### read tree
tr <- read.tree(TREE_FILE)



### short tip labels
V <- tr$tip.label

tr$tip.label <- sub('^.+RdRp(.+)$', '\\1', tr$tip.label)
tr$tip.label <- sub('^(.+)_length_.+$', '\\1', tr$tip.label)

names(V) <- tr$tip.label



### rotate tree
for (i in c(43, 52)) {
	tr <- rotate(tr, node = i)
	tr <- read.tree(text=write.tree(tr))
}



### branch colors
br_col <- rep('black', length(tr$edge))
br_col[ which.edge(tr, group = 18:42) ] <- 'royalblue'
br_col[ which.edge(tr, group = 9:17)  ] <- 'seagreen3'
br_col[ which.edge(tr, group = 2:8)   ] <- 'darkorchid3'



### tip colors
TIP_COLOR <- setNames(
				rep('black', Ntip(tr)),
				tr$tip.label
)
TIP_COLOR[ grepl('^GFD_', names(TIP_COLOR)) ] <- 'brown2'
TIP_COLOR[ 'Hubei_earwig_virus_2' ] <- 'deepskyblue'



### plot
pdf(
	'ICTV_GFD_picobirna_RdRp_tree.pdf',
	width = 18 / 2.54,
    height = 24 / 2.54
)


par(mar = rep(0.1, 4))


plot(
	tr,
	tip.color = TIP_COLOR[ tr$tip.label ],
	edge.color = br_col,
	edge.width = ifelse(br_col == 'black', 1, 3),
	x.lim = 5.5,
	font = 1
)


tree_support_pch(
	tr,
	bs_pos = 'topleft',
	bs_name = 'BP',
	max = 100
)


nodelabels(
	c('Genogroup 1', 'Genogroup 2'),
	node = c(60, 52),
	col = c('royalblue', 'seagreen3'),
	bg = NA,
	frame = 'none',
	adj = c(1.1, -1),
	cex = 1.1
)

nodelabels(
	'Genogroup 3',
	node = 45,
	col = 'darkorchid3',
	bg = NA,
	frame = 'none',
	adj = c(1.1, 2),
	cex = 1.1
)


add.scale.bar(x = 0, y = 38)


dev.off()



### subset alignment for ESPript plot
sele1 <- c('HumanPBVHy005102', grep('^GFD_', tr$tip.label[42:18], value = TRUE))
sele2 <- c('PBVhumanCDC23', grep('^GFD_', tr$tip.label[17:9], value = TRUE))
sele <- c(sele1, sele2)
SELE <- V[sele]


# read alignment
A <- read.fasta(MSA_FILE, seqtype = 'AA')


# subset alignment, rename sequences, remove empty columns
A <- A[SELE]
names(A) <- sele
A <- remove_empty_columns(A)


# write alignment
write.fasta(
	sequences = A,
	names = names(A),
	file.out = 'ICTV_GFD_picobirna_RdRp_msa_outgr.SELE.fasta'
)



### contigs order on a tree
contigs <- V[ grep('^GFD_', rev(tr$tip.label), value = TRUE) ]
contigs <- sub('_[0-9]+$', '', contigs)

write.table(
	contigs,
	row.names = FALSE,
	col.names = FALSE,
	quote = FALSE,
	file = 'picobirna_contigs_tree_order.txt'
)
