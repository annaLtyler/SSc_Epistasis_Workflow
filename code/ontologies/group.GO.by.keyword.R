#This function groups GO terms based on keywords entered by hand
#keyword.list <- list()
#keyword.list[[1]] <- c("signal transduction", "signaling")
#keyword.list[[2]] <- c("binding")
#keyword.list[[3]] <- c("development")
#keyword.list[[4]] <- c("collagen", "matrix", "cytoskeleton", "extracellular", "organization")


group.GO.by.keyword <- function(go.terms, keyword.list){
	
	term.order <- vector(mode = "list", length = length(keyword.list))
	for(i in 1:length(keyword.list)){
		all.locale <- NULL
		for(j in 1:length(keyword.list[[i]])){
			all.locale <- unique(c(all.locale, grep(keyword.list[[i]][j], go.terms)))
			}
		term.order[[i]] <- go.terms[all.locale]
		}
	names(term.order) <- names(keyword.list)

	#remove duplicate entries
	#keep the final placement of each term
	term.placement <- t(sapply(go.terms, function(y) sapply(term.order, function(x) length(which(x == y)))))
	final.placement <- apply(term.placement, 1, function(x) if(any(x == 1)){max(which(x == 1))})

	no.dups <- vector(mode = "list", length = length(keyword.list))
	names(no.dups) <- names(keyword.list)
	for(i in 1:length(final.placement)){
		if(!is.null(final.placement[[i]])){
		no.dups[[final.placement[[i]]]] <- c(no.dups[[final.placement[[i]]]], names(final.placement)[i])
		}
	}

	ordered.names <- unlist(no.dups)
	orderV <- match(ordered.names, go.terms)
	#go.terms[orderV]
	
	remaining <- setdiff(1:length(go.terms), orderV)
	if(length(remaining) > 0){
		no.dups[[length(no.dups)+1]] <- go.terms[remaining]
	}
	
	return(no.dups)
}