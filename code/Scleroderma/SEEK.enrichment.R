#This function looks at enrichment of genes correlated with
#query genes in the SEEK database
#It requires having files downloaded from http://seek.princeton.edu
#in the working directory. The files should be named for the subset
#of data sets they were gathered from.
# query.genes <- c("WNT5A", "RBMS3")


SEEK.enrichment <- function(query.genes, file.path = ".", pval = 0.05){
		
	gene.files <- list.files(path = file.path, pattern = '.txt', full.names = TRUE)
	list.types <- sapply(strsplit(gsub(".txt", "", gsub(path.expand(file.path), "", gene.files)), "/"), function(x) x[2])
	
	all.coexpression <- sig.coexpression <- thresh.coexpression <- vector(mode = "list", length = length(list.types))
	names(all.coexpression) <- names(sig.coexpression) <- names(thresh.coexpression) <- list.types
	
	all.thresholds <- rep(NA, length(list.types))
	
	enrich.sig.file <- paste0(file.path, "/Enrichment.Sig.RData")
	
	#get the functional enrichments for each gene set
	if(file.exists(enrich.sig.file)){
		all.enrich.sig <- readRDS(enrich.sig.file)
		}else{
		all.enrich.sig <- all.enrich.thresh <- vector(mode = "list", length = length(list.types))
		names(all.enrich.sig) <- names(all.enrich.thresh) <- list.types
		
		for(ty in 1:length(list.types)){
			cat(list.types[ty], "\n")
			gene.list <- read.table(list.files[ty], stringsAsFactors = FALSE, header = TRUE, sep = "\t")
			head(gene.list)
			
			all.coexpression[[ty]] <- cbind(gene.list[,"Coexpression.Score"], gene.list[,"P.Value"])
			
			sig.genes <- gene.list[which(gene.list[,"P.Value"] <= seek.pval),]
		
			sig.coexpression[[ty]] <- sig.genes
				
			if(nrow(total.genes) > 0){		
				all.enrich.sig[[ty]] <- gprofiler(c(query.genes, sig.genes[,"Gene"]), "hsapiens", hier_filtering = "strong", ordered_query = T, sort_by_structure = F, src_filter = "GO")
				}
			}
		saveRDS(all.enrich.sig, enrich.sig.file)
		}
	
	return(all.enrich.sig)
	
	
}