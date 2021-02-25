# sample_space = 3 will ensure that after each three samples there is an empty space

plot_contigs_abundance <- function (M,
                                    contig_col = NA, contig_font = NA,
                                    sample_space = NA,
                                    brackets = NA, brackets_lab = NA, brackets_name = 'RPKM') {

	### calculate y coordinates
	Y <- 1:nrow(M)
	if (!is.na(sample_space)) {

		for (i in seq_along(Y)) {

			if (i %% sample_space == 0) {

				if (i == length(Y)) { next }

				idx <- (i+1):length(Y)
				Y[idx] <- Y[idx] + 1

			}
		}
	}


	### palette
	colfunc <- colorRampPalette(c('blue', 'red'))
	pal <- data.frame(
			col = colfunc(length(brackets) - 1),
			Above = brackets[-length(brackets)],
			BelowEqual = brackets[-1],
			stringsAsFactors = FALSE
	)

	if (is.na(contig_col[1])) {
		contig_col <- rep('black', ncol(M))
	}

	if (is.na(contig_font[1])) {
		contig_font <- rep(1, ncol(M))
	}


	### plot
	plot(
		NA,
		xlim = c(0, ncol(M)),
		ylim = c(max(Y) + 3, 0),
		axes = FALSE,
		ann = FALSE
	)
	x_shift <- ncol(M) / 20
	y_shift <- nrow(M) / 50

	# sample names
	text(
		labels = rownames(M),
		x = -x_shift,
		y = Y - 0.5,
		adj = 1,
		xpd = TRUE
	)

	# contig names
	text(
		labels = colnames(M),
		col = contig_col,
        font = contig_font,
		x = 1:ncol(M) - 0.5,
		y = -y_shift,
		adj = 0,
		srt = 90,
		xpd = TRUE
	)

	# counts
	for (i in 1:nrow(M)) {

		for (j in 1:ncol(M)) {

			rect(
				xleft = j - 1,
				xright = j,
				ybottom = Y[i],
				ytop = Y[i] - 1,
				col = pal$col[pal$Above < M[i, j] & pal$BelowEqual >= M[i, j]]
			)

		}
	}

	# legend
    X <- seq(from = 0, to = ncol(M), length = length(brackets))
    rect(
        xleft = X[-length(X)],
        xright = X[-1],
        ybottom = Y[length(Y)] + 3,
        ytop = Y[length(Y)] + 2,
        col = pal$col,
        border = NA
    )

	text(
		labels = format(brackets_lab, scientific = FALSE),
		x = X[which(brackets %in% brackets_lab)],
		y = Y[length(Y)] + 3 + y_shift,
		adj = 1,
		srt = 90,
		xpd = TRUE
	)

	text(
		labels = brackets_name,
		x = -x_shift,
		y = Y[length(Y)] + 2.5,
		adj = 1,
		xpd = TRUE
	)

}
