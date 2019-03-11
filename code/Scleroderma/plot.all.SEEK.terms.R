#This function can be run after SEEK.enrichment() 
#to plot enrichment results in a single matrix
plot.all.SEEK.terms <- function(all.enrich){

	all.mats <- vector(mode = "list", length = length(all.enrich))
	names(all.mats) <- names(all.enrich)
	#get all the terms found for all gene pairs
	all.terms <- NULL
	for(gp in 1:length(all.enrich)){
		gp.enrich <- all.enrich[[gp]]
		all.terms <- unique(c(all.terms, unlist(lapply(gp.enrich, function(x) x[which(x[,"col"] != "gray"),"term.name"]))))
		}
	#put the terms in an order reflective of our original order, for functional grouping
	all.terms <- unique(all.terms[unlist(lapply(highlight.terms, function(x) grep(x, all.terms, ignore.case = TRUE)))])

	#make a matrix of all terms, and identify which were found for each gene pair and each data set
	for(gp in 1:length(all.enrich)){
		gp.enrich <- all.enrich[[gp]]	
		term.mat <- matrix(NA, nrow = length(all.terms), ncol = length(gp.enrich))
		rownames(term.mat) <- all.terms
		colnames(term.mat) <- names(gp.enrich)

		for(ds in 1:length(gp.enrich)){
			has.term <- lapply(all.terms, function(x) grep(x, gp.enrich[[ds]][,"term.name"]))
			term.p <- lapply(has.term, function(x) if(length(x) > 0){max(-log10(gp.enrich[[ds]][x,"p.value"]))}else{NA})
			term.mat[,ds] <- unlist(term.p)
			}
		if(any(!is.na(term.mat))){
			pheatmap(term.mat, main = names(all.enrich)[[gp]], cluster_cols = FALSE, cluster_rows = FALSE)
			}else{
			plot.text("No Enrichments")
			}
		all.mats[[gp]] <- term.mat
		} #end looping through gene pairs	
	return(all.mats)
	} #end function
