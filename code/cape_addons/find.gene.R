#This function takes in a gene and a data.obj
#and returns the linkage block position of
#the gene

find.gene <- function(data.obj, gene.name){
	
	all.genes <- read.table("~/Documents/Data/Mice/mouseGenes.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
	gene.locale <- which(all.genes[,3] == gene.name)
	if(length(gene.locale) == 0){
		cat("This gene could not be found.\n")
		return()
		}
	gene.chr <- all.genes[gene.locale,5]
	gene.start <- all.genes[gene.locale,6]
	gene.end <- all.genes[gene.locale,7]

	chr.locale <- which(data.obj$chromosome == gene.chr)
	marker.pos <- data.obj$marker.loc[chr.locale]
	marker.names <- data.obj$marker.names[chr.locale]
	
	
	before.start <- which(marker.pos <= gene.start)
	
	if(length(before.start) == 0){
		gene.block <- paste("Chr", gene.chr, "_", 1, sep = "")
		return(gene.block)
		}
	
	after.end <- which(marker.pos >= gene.end)
	if(length(after.end) == 0){
		block.chr <- multi.strsplit(names(data.obj$linkage.blocks.collapsed), patterns = c("Chr", "_1", "_2"))
		chr.locale <- which(block.chr == gene.chr)
		return(names(data.obj$linkage.blocks.collapsed)[max(chr.locale)])
		}
	
	
	flanking.markers <- c(max(before.start), min(after.end))
	
	marker.num <- get.marker.num(data.obj, marker.names[flanking.markers])
	
	marker.locale <- lapply(data.obj$linkage.blocks.collapsed, function(x) which(x %in% marker.num))
	block.which <- which(unlist(lapply(marker.locale, length)) > 0)
	
	return(names(block.which))
	
}