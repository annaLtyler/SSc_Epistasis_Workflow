#This function formats the output from the mouse map
#converter to go into BioMart


region2BioMart <- function(filename = "Block.Coord.txt"){
	
	
	file.base <- strsplit(filename, ".txt")[[1]][1]
	regions <- read.table(filename, stringsAsFactors = FALSE)

	#take out the covariates on chromosome 0
	regions <- regions[which(regions[,2] != 0),]
	
	region.names <- sapply(strsplit(regions[,1], "\\."), function(x) x[1])
	region.pos <- sapply(strsplit(regions[,1], "\\."), function(x) x[2])
	
	biomart.table <- NULL
	genome.browser.table <- NULL
	
	for(bl in 1:length(region.names)){
		#find the min and max for each block
		block.min.locale <- intersect(which(region.names == region.names[bl]), which(region.pos == "min"))
		block.max.locale <- intersect(which(region.names == region.names[bl]), which(region.pos == "max"))
		block.chr <- regions[block.min.locale,2]
		block.begin <- as.numeric(regions[block.min.locale,3])
		block.end <- as.numeric(regions[block.max.locale,3])
		#add both plus and minus strand, since we only know the
		#locations of our genes.
		table.row <- c(block.chr, block.begin, block.end)
		biomart.table <- rbind(biomart.table, table.row)
		genome.browser.table <- rbind(genome.browser.table, paste(block.chr, ":", block.begin, "..", block.end, collapse = "", sep = ""))
		}
	
	write.table(biomart.table, paste(file.base, ".BioMart.txt", sep = ""), quote = FALSE, sep = ":", col.names = FALSE, row.names = FALSE)
	write.table(genome.browser.table, paste(file.base, ".GenomeBrowser.txt", sep = ""), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
	
}