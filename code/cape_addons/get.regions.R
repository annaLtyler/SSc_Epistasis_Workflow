#This function finds regions surounding marker blocks


get.regions <- function(data.obj, p.or.q = 0.05, r2.thresh = 0.8, region.buffer = 1){
	
	
	old.r2.thresh <- data.obj$r2.thresh
	if(is.null(old.r2.thresh) || old.r2.thresh != r2.thresh){
		cat("running linkage.blocks() to determine block structure...\n")
		data.obj <- linkage.blocks(data.obj, p.or.q = p.or.q, r2.thresh = r2.thresh)
		}

	linkage.blocks <- data.obj$linkage.blocks.collapsed

	get.region <- function(block){
		marker.loc <- which(colnames(data.obj$geno) %in% block)
		marker.ch <- unique(data.obj$chromosome[marker.loc])
		marker.pos <- as.numeric(data.obj$marker.location[marker.loc])
		min.pos <- max(c(min(marker.pos)-region.buffer, 0))
		max.pos <- max(marker.pos) + region.buffer
		result <- c(marker.ch, min.pos, max.pos)
		names(result) <- c("chromosome", "min.pos", "max.pos")
		return(result)
		}
	
	block.regions <- t(sapply(linkage.blocks, get.region))
	block.regions <- apply(block.regions, 2, as.numeric)
	rownames(block.regions) <- names(linkage.blocks)
	
	data.obj$block.regions <- block.regions
	return(data.obj)
}