#This function tabulates interactions and main effects by phenotype

tabulate.results <- function(data.obj, specific.markers = NULL){
	
	collapsed.net <- data.obj$collapsed.net
	
	pheno.names <- colnames(data.obj$pheno)
	pheno.locale <- match(pheno.names, colnames(collapsed.net))
	
	all.col <- 1:dim(collapsed.net)[2]
	if(length(pheno.locale) > 0){
		just.int.locale <- all.col[-pheno.locale]
		}
	
	just.int <- collapsed.net[just.int.locale, just.int.locale]
	just.pheno <- collapsed.net[,pheno.locale]
	
	int.table <- matrix(NA, nrow = 3, ncol = 1)
	colnames(int.table) <- "interactions"
	int.table[1,1] <- length(which(just.int > 0))
	int.table[2,1] <- length(which(just.int < 0))
	int.table[3,1] <- neg.int + pos.int
	
	pheno.table <- matrix(NA, nrow = 3, ncol = length(pheno.locale))
	rownames(pheno.table) <- c("positive", "negative", "total")
	colnames(pheno.table) <- pheno.names

	for(i in 1:length(pheno.locale)){
		pheno.table[1,i] <- length(which(just.pheno[,i] > 0))
		pheno.table[2,i] <- length(which(just.pheno[,i] < 0))
		pheno.table[3,i] <- length(which(just.pheno[,i] != 0))
		}

	marker.table <- matrix(NA, nrow = 3, ncol = (length(specific.markers)*3))
	header.names <- NULL
	for(i in 1:length(specific.markers)){
		header.names <- c(header.names, paste(specific.markers[i], c("-incoming", "-outgoing", "-total"), sep = ""))
		}
	colnames(marker.table) <- header.names
	rownames(marker.table) <- rownames(pheno.table)
	
	if(!is.null(specific.markers)){
		for(i in 1:length(specific.markers)){
			marker.col <- which(colnames(just.int) == specific.markers[i])
			marker.row <- which(rownames(just.int) == specific.markers[i])

				incoming <- just.int[,marker.col]
				marker.table[1,((i*3)-2)] <- length(which(incoming > 0))	
				marker.table[2,((i*3)-2)] <- length(which(incoming < 0))
				marker.table[3,((i*3)-2)] <- length(which(incoming != 0))


				outgoing <- just.int[marker.row,]
				marker.table[1,((i*3)-1)] <- length(which(outgoing > 0))	
				marker.table[2,((i*3)-1)] <- length(which(outgoing < 0))
				marker.table[3,((i*3)-1)] <- length(which(outgoing != 0))				


				all.int <- c(just.int[marker.row,], just.int[,marker.col])
				marker.table[1,(i*3)] <- length(which(all.int > 0))
				marker.table[2,(i*3)] <- length(which(all.int < 0))
				marker.table[3,(i*3)] <- length(which(all.int != 0))
			}
		}else{
			marker.table <- NULL
		}

	final.table <- cbind(int.table, pheno.table, marker.table)
	
	write.table(final.table, "Effects.Counts.txt", sep = "\t", quote = FALSE)
	
	
	
	
}