#This function uses significant SNPs to define linkage blocks


all.linkage.blocks <- function(data.obj, r2.thresh = 0.8){
	

	marker.geno <- data.obj$geno
	sorted.markers <- colnames(marker.geno)
	
	data.obj$r2.thresh <- r2.thresh

	all.cor <- cor(marker.geno, use = "pairwise.complete.obs")^2
	#zero out the diagonal and the lower triangle
	all.cor[lower.tri(all.cor, diag = TRUE)] <- 0
	
	in.ld <- which(all.cor >= r2.thresh, arr.ind = TRUE)


	#start with a linkage block list that contains one marker per block
	linkage.blocks <- as.vector(sorted.markers, mode = "list")
	names(linkage.blocks) <- paste("Block", 1:length(linkage.blocks), sep = "")
	
	#if there are no markers in LD, just return the list as is.
	if(length(in.ld) == 0){
			data.obj$linkage.blocks <- linkage.blocks
			return(data.obj)
			}
	
	#otherwise, go through the in.ld matrix and combine markers that are linked
	linked.markers <- list(in.ld[1,])
	
	if(dim(in.ld)[1] > 1){
		for(i in 2:length(in.ld[,1])){

			row.to.check <- in.ld[i,]
			#look for blocks that already contain markers we're looking at
			#The lapply function subtracts the length of the combined unique
			#vector from the length of the total vector. If the result is 
			#positive, there are common nodes between the new row and an
			#existing block
			common.nodes <- lapply(linked.markers, function(x) length(c(x,row.to.check))-length(unique(c(x,row.to.check))))
			shared.node.locale <- which(common.nodes > 0)

			#if there are no shared nodes, start a new block 
			if(length(shared.node.locale) == 0){
				linked.markers[[(length(linked.markers)+1)]] <- as.vector(row.to.check)
				}

			#if we need to combine the new row with one other block
			if(length(shared.node.locale) == 1){
				linked.markers[[shared.node.locale]] <- unique(c(linked.markers[[shared.node.locale]], row.to.check))
				}
		
			if(length(shared.node.locale) > 1){
				message("There is more than one common node")
				}

			}
		}

	#go through the linkage blocks and adjust the marker block list
	for(i in 1:length(linked.markers)){
		marker.names <- colnames(marker.geno)[sort(linked.markers[[i]])]
		# cat(i, marker.names, "\n")
		list.locale <- match(marker.names, linkage.blocks)
		#replace the first marker in the block with all markers
		linkage.blocks[[min(list.locale)]] <- marker.names
		#nullify the other blocks
		to.nullify <- list.locale[-1]
		for(j in 1:length(to.nullify)){
			# cat("\ttaking out", linkage.blocks[[to.nullify[j]]], "\n")
			linkage.blocks[[to.nullify[j]]] <- 0
			}
		}
	#remove all the blocks that have been converted to 0
	to.nullify <- which(linkage.blocks %in% 0)
	while(length(to.nullify) > 0){
		linkage.blocks[[to.nullify[1]]] <- NULL
		to.nullify <- which(linkage.blocks %in% 0)
		}
	
	
	names(linkage.blocks) <- paste("Block", 1:length(linkage.blocks), sep = "")
	

	data.obj$linkage.blocks <- linkage.blocks
		
	return(data.obj)
	
}