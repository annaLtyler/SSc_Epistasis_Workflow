#This script reads in Sanger SNP database output
#and compares the SNPs to a gene list. If gene.list 
#is NULL, it returns all the genes with SNPs
#the strain is the strain that is not the reference
#strain (B6) in which you want to find SNPs

genes.with.snps <- function(filename, gene.list = NULL, strain = "C3HHeJ"){
	
	base.file <- strsplit(filename, ".csv")[[1]][1]
	all.snps <- as.matrix(read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE))
	strain.col <- which(colnames(all.snps) == strain)

	snp.locale <- which(all.snps[,strain.col] != ".")
	u_genes <- unique(all.snps[snp.locale,"Gene"])
	
	if(is.null(gene.list)){
		cat(length(u_genes), "genes with SNPs in", filename)
		write.table(u_genes, paste(base.file, ".with.SNPs.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
		}else{
		overlap <- intersect(gene.list, u_genes)
		cat(length(overlap), "genes in gene.list with SNPs in", filename)
		write.table(overlap, paste(base.file, ".with.SNPs.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
		}
	
}