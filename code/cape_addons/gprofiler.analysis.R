#This script takes in a list of genes, does a gprofiler 
#analysis and writes out the table filtered by a
#specified (T+Q)/T
#It writes out the file of significant terms
#and returns the associated go.terms for use in
#MouseMine

gprofile.analysis <- function(gene.list, recall.thresh = 0.5, filename = "gprofiler.list.txt"){
	
	library(gProfileR)
	gene.summ <- gprofiler(gene.list, organism = "mmusculus")
	high.recall.locale <- which(gene.summ[,"recall"] >= 0.5)
	high.recall <- gene.summ[high.recall.locale,]
	write.table(high.recall, filename, quote = FALSE, sep = "\t", row.names = FALSE)
	
	go.terms <- unique(high.recall[,"term.id"])
	return(go.terms)
	
	
	
}