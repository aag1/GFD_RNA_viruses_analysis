### A - alignment read by 'seqinr' package: 
### a list where each element is a sequence represented by vector of characters


remove_empty_columns <- function (A) {

	bad <- sapply(1:length(A[[1]]), function (i) {

				column <- unlist(lapply(A, function (v) v[i]))
				all(column == '-')

	})

	A <- lapply(A, function (v) v[!bad])

	return(A)

}
