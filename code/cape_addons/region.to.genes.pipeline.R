#This is a pipeline that goes from a single Source-Target
#relationship to candidate genes

#========================================================================================
#load the data object you want to examine
#========================================================================================

region.to.genes.pipeline <- function(data.obj, num.cand.interactions = 30, r2.thresh = 0.8){
	
	top.interactor.info(data.obj, num.pairs = num.cand.interactions, r2.thresh = r2.thresh)
	cat(paste("There is new output: Top.", num.cand.interactions, ".Interactions.r2.thresh.", r2.thresh, ".txt", sep = ""), "\n")	
	
	cat("\nUse this file to select one Source block and one Target block to focus on\n")
	source.block <- readline("Using Interaction.Plots.pdf, select one block to be the Source block:\n")
	target.block <- readline("Using Interaction.Plots.pdf, select one block to be the Target block:\n")	

	
	cat("\nTo translate marker regions to bp we need the bp location of each marker.\n")
	cat("This file should be named marker_bp.txt.\n")
	all.set <- readline(prompt = "Press enter when ready (x to quit):\n")
	
	if(all.set == "x"){
		return()
		}
	

	write.block.coord(data.obj, r2.thresh = r2.thresh, block.num = as.numeric(source.block), output.file = "Block.Coord.Source.txt")
	write.block.coord(data.obj, r2.thresh = r2.thresh, block.num = as.numeric(target.block), output.file = "Block.Coord.Target.txt")	
	
	# cat("Convert the Block.Coord files to bp using http://cgd.jax.org/mousemapconverter/\n")
	# cat("Name the result files\n\tBlock.Coord.bp.Source.txt\n\tBlock.Coord.bp.Target.txt\n")
	# all.set <- readline(prompt = "Press return when ready.\n")
	
	region2BioMart(filename = "Block.Coord.Source.txt")
	region2BioMart(filename = "Block.Coord.Target.txt")
	
	cat("\n\nUse These files in BioMart\n\tBlock.Coord.bp.Source.BioMart.txt\n\tBlock.Coord.bp.Target.BioMart.txt\n")
	cat("\n\nUse These files in the Mouse Genome Browser\n\tBlock.Coord.bp.Source.GenomeBrowser.txt\n\tBlock.Coord.bp.Target.GenomeBrowser.txt\n")
	
	cat("\nBioMart is here:\n\thttp://biomart.informatics.jax.org/")
	cat("\n\nThe genome browser is here:\n\thttp://gbrowse.informatics.jax.org/cgi-bin/gb2/gbrowse/mousebuild38/")
	
	cat("\n\nYou can look up SNPs in the regions by going here:\n\thttp://www.sanger.ac.uk/sanger/Mouse_SnpViewer/rel-1303\n\tThe format for locations is  chr:startbp-endbp.\n")
	cat("Rename output files to:\n\tSNPs_source.csv\n\tSNPs_target.csv")
	
	all.set <- readline("Press Enter when files are ready (x to quit, or n to skip to gold standard comparison):\n")
	
	if(all.set == "x"){
		return()
		}
	
	if(all.set != "n"){
	
		genes.with.snps(filename = "SNPs_source.csv")
		cat("\n")
		genes.with.snps(filename = "SNPs_target.csv")
		cat("\n")
		}else{
		source.biomart.genes <- read.table("mart_export_source.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
		source.genes <- unique(source.biomart.genes[,"Associated.Gene.Name"])
		write.table(source.genes, "SNPs_source.with.SNPs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		target.biomart.genes <- read.table("mart_export_target.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
		target.genes <- unique(target.biomart.genes[,"Associated.Gene.Name"])	
		write.table(target.genes, "SNPs_target.with.SNPs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		}
	
	
		
	cat("\nYou can now compare these gene lists with gold standard genes.\n")
	cat("Make sure gold standard files have the file name structure gold_standard_phenotype.txt\n")
	
	all.set <- readline("To continue and compare the gene lists to a gold standard, press enter (x to quit)\n")
	
	if(all.set == "x"){
		return()
		}
	
	
	gold.files <- list.files(pattern = "gold_standard")
	if(length(gold.files) == 0){
		stop("Please make sure gold standard files are labeled with 'gold_standard'\n")
		}
		
	cat("I will compare the source and target genes to the following gold standard lists:\n")
	cat(gold.files, sep = "\n")

	source.genes <- as.matrix(read.table("SNPs_source.with.SNPs.txt", stringsAsFactors = FALSE))
	target.genes <- as.matrix(read.table("SNPs_target.with.SNPs.txt", stringsAsFactors = FALSE))

	for(i in 1:length(gold.files)){
		compare.to.gold(gold.files[i], source.genes, "Source")
		compare.to.gold(gold.files[i], target.genes, "Target")
		}

	
	cat("If there are gold standard genes in the source and target regions, try using IMP to see if they interact.\n")
	cat("Otherwise, you may want to use IMP or PILGRM to expand your candidate lists.\n")
	
	
	

	
}