#This function uses the "is_a" designation to get
#parents of a given term
#obo.list <- readRDS("~/Documents/Data/Ontologies/GO.RData")


get.parents <- function(term.name, obo.list, out.type = c("name", "ID"), unstructured = TRUE){

	if(term.name == "not found" || is.na(term.name) || length(term.name) == 0){return(NA)}

	obo.name <- sapply(obo.list, function(x) x[which(x[,1] == "name:")[1],2])
	id <- sapply(obo.list, function(x) x[1,2])
	isa <- lapply(obo.list, function(x) x[which(x[,1] == "is_a:"),2])
	term.names <- cbind(obo.name, id)

	term.locale <- unlist(apply(term.names, 2, function(x) which(x == term.name)))

	if(length(term.locale) == 0){return(NA)}
	
	if(out.type == "name"){	
		term.isa <- lapply(strsplit(isa[[term.locale]], " ! "), function(x) x[2])
		}else{
		term.isa <- lapply(strsplit(isa[[term.locale]], " ! "), function(x) x[1])
		}
	
	one.up <- function(term.locale, type){
		if(length(term.locale) == 0){return(NULL)}
		next.up <- isa[[term.locale]]
		if(length(next.up) == 0){
			return(NULL)
		}else{
			if(type == "name"){
				isa.term <- lapply(next.up, function(x) strsplit(x, " ! ")[[1]][2])
			}else{
				isa.term <- lapply(next.up, function(x) strsplit(x, " ! ")[[1]][1])
			}
		return(isa.term)	
		}
	}

	
	idx <- 1
	all.terms <- list()
	all.terms[[idx]] <- term.isa
	while(!all(is.null(unlist(term.isa)))){
		idx <- idx + 1
		term.locale <- lapply(term.isa, function(x) unlist(apply(term.names, 2, function(y) which(y == x))))
		term.isa <- lapply(term.locale, function(x) one.up(x, out.type))[[1]]
		if(!all(is.null(unlist(term.isa)))){
			all.terms[[idx]] <- term.isa
		}
	}

if(unstructured){
	return(as.vector(unlist(all.terms)))
}else{
	return(all.terms)
}
}
