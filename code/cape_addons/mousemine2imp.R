#This function takes in a file from a mousemine query and
#gets it ready to put into IMP. It filters the gene list
#down to unique genes and puts the pubmed IDs in a single row

mousemine2imp <- function(filename, max.ids = NULL){
	
	file.base <- strsplit(filename, ".txt")[[1]][1]
	
	genes <- as.matrix(read.table(filename, sep = "\t", stringsAsFactors = FALSE))
	u_genes <- unique(genes[,2])
	
	gene.list <- vector(mode = "list", length = length(u_genes))
	names(gene.list) <- u_genes
	
	for(i in 1:length(u_genes)){
		gene.locale <- which(genes[,2] == u_genes[i])
		pubmed.ids <- genes[gene.locale,5]
		gene.list[[i]] <- unique(pubmed.ids[which(!is.na(pubmed.ids))])
		}
	
	if(is.null(max.ids)){
		max.len <- max(sapply(gene.list, length))
		}else{
		max.len <- max.ids	
		}
	
	imp.table <- matrix(NA, nrow = length(u_genes), ncol = max.len)
	rownames(imp.table) <- u_genes
	
	
	if(is.null(max.ids)){
		for(i in 1:length(gene.list)){
			imp.table[i,1:length(gene.list[[i]])] <- gene.list[[i]]
			}
		}else{
		for(i in 1:length(gene.list)){
			imp.table[i,] <- gene.list[[i]][1:max.len]
			}
		}
	
	write.table(imp.table, file = paste(file.base, ".imp.txt", sep = ""), quote = FALSE, sep = "\t", col.names = FALSE, na = "")
	
}