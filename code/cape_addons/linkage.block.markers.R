#This function writes out the markers in each linkage block


linkage.block.markers <- function(data.obj, r2.thresh = 0.8){
	
	orig.r2 <- cross$r2.thresh
	
	if(is.null(orig.r2) || orig.r2 != r2.thresh){
		data.obj <- all.linkage.blocks(data.obj, r2.thresh)
		}
	
	blocks <- data.obj$linkage.blocks
	
	get.markers <- function(block){
		marker.locale <- which(colnames(data.obj$geno) %in% block)
		marker.names <- data.obj$marker.names[marker.locale]
		return(marker.names)
		}
	
	marker.blocks <- lapply(blocks, get.markers)
	
	max.len <- max(sapply(marker.blocks, length))
	marker.mat <- matrix(NA, nrow = length(blocks), ncol = max.len)
	
	for(i in 1:length(marker.blocks)){
		marker.mat[i,1:length(marker.blocks[[i]])] <- marker.blocks[[i]]
		}
	
	rownames(marker.mat) <- names(blocks)
	
	write.table(marker.mat, paste("Linkage.Block.Membership.r2.thresh.", r2.thresh, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = FALSE, na = "")
	
	invisible(data.obj)
}