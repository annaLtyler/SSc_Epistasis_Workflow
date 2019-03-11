#This function looks through a snp file and a file downloaded
#from BioMart listing the genes in given regions
#The function tries to filter the gene list based on which
#have SNPs in them.
#gene.file <- "mart_export-1.csv"; snp.file <- "block1_snps.csv"

filter.genes.by.SNP <- function(snp.file, gene.file){

	snps <- as.matrix(read.table(snp.file, sep = ",", stringsAsFactors = FALSE, header = TRUE))
	genes <- as.matrix(read.table(gene.file, sep = ",", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)	)
	
	
	get.gene.snps <- function(gene.row){
		
		gene.chr <- as.numeric(gene.row["Chromosome.Name"])
		start.pos <- as.numeric(gene.row["Gene.Start..bp."])
		end.pos <- as.numeric(gene.row["Gene.End..bp."])
		
		all.chr.locale <- which(snps[,"Chr"] == gene.chr)
		snp.table <- snps[all.chr.locale,]
		snp.locale <- intersect(which(as.numeric(snp.table[,"Position"]) >= start.pos), which(as.numeric(snp.table[,"Position"]) <= end.pos))
		if(length(snp.locale) > 1){
			new.gene <- c(gene.row[c("VEGA.gene.ID", "Description", "Chromosome.Name", "Gene.Start..bp.", "Gene.End..bp.", "Strand", "Gene.Biotype", "Status..gene.", "Status..transcript.", "Transcript.Biotype", "External.Gene.ID")], length(snp.locale))
			names(new.gene) <- NULL
			return(new.gene)
			}
		}
	
	
	genes.with.snps <- apply(genes, 1, get.gene.snps)
	snp.table <- matrix(unlist(genes.with.snps), ncol = 12, byrow = TRUE)
	
	write.table(snp.table, file = paste("Genes.with.SNPS.", snp.file, sep = ""), quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
	
	
	
}