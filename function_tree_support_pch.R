# max - should be either 1 or 100, indicating that support values are in [0,1] or [0,100] range, correspondingly

tree_support_pch <- function (tr, bs_pos, bs_name, max) {

	n <- which(!(is.na(tr$node.label)) & tr$node.label!='NA' & tr$node.label!='')

	v <- tr$node.label[n]
	v <- as.numeric(v)

	if (max==1) { v <- v*100 }

	bg <- sapply(v, function (x) {
					if (x>=90) {return("black")}
					if (x<90 && x>=70) {return("grey75")}
					if (x<70) {return("white")}
	})

	nodelabels(node=n+Ntip(tr), pch=21, bg=bg, col="black")

	if (!is.na(bs_pos)) {
		legend( bs_pos,
			legend=eval(substitute(expression("90%" <= bs, "70%" <= bs * " < 90%", bs < "70%"), list(bs=bs_name))),
			pch=21, pt.bg=c("black","grey","white"),
			cex=1, xpd=TRUE )
	}

}
