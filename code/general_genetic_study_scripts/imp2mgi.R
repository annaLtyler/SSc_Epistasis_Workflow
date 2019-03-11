#This script makes a query for MGI out of the IMP network
#downloaded weights


imp2mgi <- function(filename, threshold = 0.9){
	impnet <- read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)
	above.thresh <- which(impnet[,3] >= threshold)
	ugenes <- unique(c(impnet[above.thresh,1], impnet[above.thresh,2]))
	write.table(ugenes, paste("MGI_", threshold, "_", filename, sep = ""), quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)	
}