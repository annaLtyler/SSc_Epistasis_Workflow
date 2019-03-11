#This function pulls out the unique genes from the candidate regions

candidate.genes <- function(filename, output.file = "just.genes.txt"){
	
	data.set <- as.matrix(read.table(filename, sep = "\t", fill = TRUE, stringsAsFactors = FALSE))
	u_genes <- unique(data.set[,2])
	write.table(u_genes, output.file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
	
}