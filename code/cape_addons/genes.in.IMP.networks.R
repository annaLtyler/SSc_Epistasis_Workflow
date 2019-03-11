#This script takes in a list of filenames of IMP networks
#and returns the list of unique genes.

genes.in.IMP.networks <- function(filenames){
	
	all.genes <- NULL
	
	for(i in 1:length(filenames)){
		edges <- as.matrix(read.table(filenames[i], sep = ",", header = TRUE, stringsAsFactors = FALSE))
		all.genes <- unique(c(all.genes, edges[,1], edges[,2]))
		}
	
	
	return(all.genes)
	
}