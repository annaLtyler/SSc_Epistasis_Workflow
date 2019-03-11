#This function compares a supplied gene list with
#genes in a gold standard output from mousemine

compare.to.gold <- function(gold.file, compare.list, description = "Source"){
	
	
		base.file <- strsplit(gold.file, ".txt")[[1]][1]
		pheno.name <- strsplit(base.file, "gold_standard_")[[1]][2]

		gold.info <- as.matrix(read.table(gold.file, stringsAsFactors = FALSE, sep = "\t", fill = TRUE))
		gold.genes <- unique(gold.info[,2])
		overlap <- intersect(compare.list, gold.genes)

		cat(length(overlap), gold.file, " gene(s) in the source region.\n")
		if(length(overlap) > 0){
			write.table(overlap, paste("Overlap_", pheno.name, "_", description, ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
			}

	
}