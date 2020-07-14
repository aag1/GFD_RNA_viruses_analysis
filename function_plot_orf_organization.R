plot_orf_organization <- function (tab, xmax = NA, scale = NA) {

    ### canvas

    if (is.na(xmax)) { xmax <- tab[1, 2] }

    par(mar = rep(0, 4))
    
    plot(
        NA,
        xlim = c(-0.1, 1.1) * xmax,
        ylim = c(-0.5, 4.5),
        xaxs = 'i', yaxs = 'i', 
        axes = FALSE, ann = FALSE, bty = 'n'
    )
    


 	### ORFs
	for (i in 2:nrow(tab)) {

		from <- tab[i, 1]
		to <- tab[i, 2]
        name <- tab[i, 3]

		utr <- from - tab[2, 1]               # ORFs defined relative to the 5'-terminal ORF
		if (utr %% 3 == 0)       { Y <- 1 }   #  0 frame
		if ((utr - 1) %% 3 == 0) { Y <- 0 }   # +1 frame
		if ((utr + 1) %% 3 == 0) { Y <- 2 }   # -1 frame

		rect(
            xleft = from,
            xright = to,
            ybottom = Y,
            ytop = Y + 1,
            col = 'grey90'
        )

        text(
            x = (from + to) / 2,
            y = Y + 0.5,
            labels = name
        )

    }



	### genome
	rect(
        xleft = tab[1, 1],
        xright = tab[1, 2], 
        ybottom = 0,
        ytop=3,
        col=NA
    )

	text(
        x = -0.01 * tab[1, 2],
        y = c(0.5, 1.5, 2.5),
        labels = c('+1', ' 0', '-1'),
        adj = 1
    ) 



    ### scale
    if (!is.na(scale)) {

        axis(
            side = 3,
            pos = 3.5,
            at = tab[1, 2] - c(scale, 0),
            labels = FALSE
        )

	   text(
            x = tab[1, 2] - c(scale, 0),
            y = 4,
            labels = c('0', paste(scale, 'nt'))
       )

    }
}
