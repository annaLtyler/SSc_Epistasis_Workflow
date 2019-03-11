#This function matches up known and predicted protein-protein interactions
#with cape interactions
#It requires a file output from the String database. The interaction
#confidence should be filtered prior to downloading. This function 
#takes everything.
# string.filename <- "~/Documents/Data/Little_Cross/little_cross_data_cape/Body_Comp_IGF_New_P/igf1.string.interactions.txt"
#When you download the file be sure to uncomment the header row

string.block.interactions <- function(data.obj, gene.name, string.filename, gene.block.name = NULL, confidence = 0.9){
	
	library(gProfileR)

	gene.interactions <- read.table(string.filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	conf.locale <- which(gene.interactions[,15] >= confidence)
	gene.interactions <- gene.interactions[conf.locale,,drop=FALSE]
	
		
	gene.locale <- which(gene.interactions[,1:2] == gene.name, arr.ind = TRUE)
	
	if(length(gene.locale) == 0){
		stop("I couldn't find the gene. Please check spelling and capitalization.")
	}
	just.gene.int <- gene.interactions[gene.locale[,1],]
	
	u_genes <- unique(c(just.gene.int[,1], just.gene.int[,2]))
	u_genes <- u_genes[-which(u_genes == gene.name)]
	
	registerDoParallel()
	gene.blocks = foreach(g = u_genes) %dopar% {
		find.gene(data.obj, g)
		}
	
	if(is.null(gene.block.name)){
		gene.block.name <- find.gene(data.obj, gene.name)
		}
	
	collapsed.net <- data.obj$collapsed.net
	int.dim <- min(dim(collapsed.net))
	block.locale <- which(rownames(collapsed.net) == gene.block.name)
	
	just.int <- collapsed.net[1:int.dim, 1:int.dim]
	
	block.targets <- names(which(just.int[block.locale,] != 0))
	targets.block <- names(which(just.int[,block.locale] != 0))
	all.int <- unique(c(block.targets, targets.block))
	
	protein.interactors <- vector(mode = "list", length = length(all.int))
	names(protein.interactors) <- all.int
	
	#find the genes in the cape interactors
	for(i in 1:length(all.int)){
		gene.block.locale <- which(unlist(lapply(gene.blocks, function(x) length(which(x == all.int[i])))) > 0)
		block.genes <- u_genes[gene.block.locale]
		if(length(block.genes) == 0){block.genes <- "none"}
		protein.interactors[[i]] <- block.genes
		}

	results.filename <- paste("Genes.In.Blocks.Interacting.With.", gene.name, ".txt", sep = "")
	int.table <- list2Matrix(protein.interactors)
	write.table(int.table, results.filename, quote = FALSE, col.names = FALSE, na = "", sep = "\t")
	
	enrichment <- vector(mode = "list", length = (length(protein.interactors)+1))
	names(enrichment) <- c(names(protein.interactors), "all")
	for(i in 1:length(protein.interactors)){
		if(protein.interactors[[i]][1] != "none"){
			enrichment[[i]] <- gprofiler(c(gene.name, protein.interactors[[i]]), organism = "mmusculus")
			}else{
			enrichment[[i]] <- "none"
			}
		}
	all.genes <- unlist(protein.interactors)
	all.genes <- all.genes[which(all.genes != "none")]
	enrichment[[length(enrichment)]] <- gprofiler(c(gene.name, all.genes), organism = "mmusculus")
	
	
	with.entries <- as.vector(which(enrichment != "none"))
	enrichment.mat <- list2Matrix(enrichment[with.entries], preserve.dim = TRUE)
	write.table(enrichment.mat, paste("Genes.In.Blocks.Interacting.With.", gene.name, ".Enrichment.txt", sep = ""), quote = FALSE, sep = "\t")

}