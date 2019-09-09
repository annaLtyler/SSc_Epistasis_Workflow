#This function creates a DAG out of the MPO
# onto.list <- readRDS("/Users/atyler/Documents/Ontologies/MPO.data.RData")

get.term.path <- function(onto.list, term.name = NULL, term.id = NULL){
	
	all.term.names <- unlist(lapply(onto.list, function(x) attr(x, "name:")))
	all.term.ids <- unlist(onto.list)
	
	if(!is.null(term.name)){
		term.locale <- which(all.term.names == term.name)
		full.path <- term.name
		}else{
		term.locale <- which(all.term.ids == term.id)	
		full.path <- term.id
		}

	if(length(term.locale) == 0){
		return(full.path)
		}
		
	current.term <- onto.list[[term.locale]]
	parent.term <- attr(current.term, "is_a")
	
	while(!is.null(parent.term)){
		new.term.split <- strsplit(parent.term, " ! ")
		new.term.id <- new.term.split[[1]][1]
		new.term.name <- new.term.split[[1]][2]
		if(!is.null(term.name)){
			full.path <- c(full.path, new.term.name)
			new.term.locale <- which(all.term.names == new.term.name)
			}else{
			full.path <- c(full.path, new.term.id)
			new.term.locale <- which(all.term.ids == new.term.id)
			}
		current.term <- onto.list[[new.term.locale]]
		parent.term <- attr(current.term, "is_a")
		}
	
	return(full.path)
}