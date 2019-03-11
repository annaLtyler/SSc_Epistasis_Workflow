#This function writes out the marker names in each linkage block

write.linkage.blocks <- function(data.obj, filename = "Linkage.Blocks.txt", pos.filename = "Linkage.Blocks.Positions.txt", print.internal.markers = TRUE, sig.figs = 3){
	
	blocks <- data.obj$linkage.blocks.collapsed
	
	if(is.null(blocks)){
		stop("linkage.blocks() or get.network() must be run before writing the linkage blocks.\n")
		}
		
	marker.names <- lapply(blocks, function(x) get.marker.name(data.obj, x))
	marker.pos <- lapply(blocks, function(x) get.marker.location(data.obj, x))
	
	if(!print.internal.markers){
		#take out all markers except first and last of each block
		stripped.marker.names <- lapply(marker.names, function(x) c(x[1], x[length(x)]))
		stripped.marker.pos <- lapply(marker.pos, function(x) c(x[1], x[length(x)]))
		}else{
		stripped.marker.names <- marker.names
		stripped.marker.pos <- marker.pos	
		}
	
	max.len <- max(sapply(stripped.marker.names, length))
	
	name.mat <- matrix(NA, ncol = max.len, nrow = length(blocks))
	rownames(name.mat) <- names(blocks)

	pos.mat <- matrix(NA, ncol = max.len, nrow = length(blocks))
	rownames(pos.mat) <- names(blocks)

	
	for(i in 1:length(blocks)){
		name.mat[i,1:length(stripped.marker.names[[i]])] <- stripped.marker.names[[i]]
		pos.mat[i,1:length(stripped.marker.pos[[i]])] <- signif(stripped.marker.pos[[i]], sig.figs)
		}
	
	
	
	write.table(name.mat, file = filename, sep = "\t", quote = FALSE, na = "", col.names = FALSE)
	write.table(pos.mat, file = pos.filename, sep = "\t", quote = FALSE, na = "", col.names = FALSE)
}