#This function writes out the marker names and positions from a data obj

marker.name.pos <- function(data.obj){
	
	
	marker.names <- data.obj$marker.names
	marker.chr <- data.obj$chromosome
	marker.pos <- data.obj$marker.location

	
	marker.table <- cbind(marker.names, marker.chr, marker.pos)
	write.table(marker.table, "Marker.Names.and.Positions.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
	
	
}