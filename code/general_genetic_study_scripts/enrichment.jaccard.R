enrichment.jaccard <- function(enrichment){
	
	overlap <- enrichment[,"overlap.size"]
	nterm <- enrichment[,"term.size"]
	nquery <- enrichment[,"query.size"]
	
	jacc <- overlap/(nterm + nquery - overlap)
	return(jacc)
	}