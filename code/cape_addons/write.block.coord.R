#This function finds genomic locations for all significant markers
#it returns single locations for individual markers and the bounds
#for linkage blocks
#the regions can be expanded in MouseMine


write.block.coord <- function(data.obj, r2.thresh = 0.8, marker.pos.file = "marker_bp.txt", block.num = NULL, output.file = "Block.Coord.txt", expand.by.bp = 0){

	#This is the maximum length of each mouse chromosome
	#obtained from http://www.informatics.jax.org/mgihome/other/mouse_facts1.shtml
	#on May 1, 2013.
	mouse.chr <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
	max.mbp <- c(197, 182, 160, 156, 153, 150, 153, 132, 124, 130, 122, 121, 120, 125, 103, 98, 95, 91, 61)
	
	
	if(is.character(block.num)){
		stop("block.num must be numeric.")
		}
	
	
	marker.bp <- read.table(marker.pos.file, stringsAsFactors = FALSE, sep = "\t")
	#make sure marker.bp is in the same order as the marker names
	marker.order <- unlist(apply(matrix(data.obj$marker.names, ncol = 1), 1, function(x) which(marker.bp[,1] == x)))
	marker.bp <- marker.bp[marker.order,]
	
	
	#get the linkage blocks based on the significant markers
	data.obj <- all.linkage.blocks(data.obj, r2.thresh = r2.thresh)
	
	blocks <- data.obj$linkage.blocks
	
	#get the genome coordinates for each block by
	#popping out to the next adjacent marker
	
	get.coord <- function(block){
		# print(block)
		marker.locale <- which(colnames(data.obj$geno) %in% block)
		marker.ch <- unique(data.obj$chromosome[marker.locale])
		if(marker.ch == 0){
			result <- c(0, 0, 0)
			names(result) <- c("chromosome", "min.boundary", "max.boundary")
			return(result)
			}
		min.locale <- min(marker.locale)
		max.locale <- max(marker.locale)
		min.pos <- marker.bp[min.locale,2]
		max.pos <- marker.bp[max.locale,2]

		min.pos <- as.numeric(min.pos) - expand.by.bp
		max.pos <- as.numeric(max.pos) + expand.by.bp
		
		# up.marker <- min.locale - 1
		# down.marker <- max.locale + 1
		# up.marker.ch <- data.obj$chromosome[up.marker]
		# down.marker.ch <- data.obj$chromosome[down.marker]

		# max.chr.len <- max.mbp[which(mouse.chr == marker.ch)]
		# if(length(max.chr.len) == 0){
			# max.chr.len <- 0
			# }

		# #if the upstream marker is on the same chromosome as the markers in question...
		# if(!is.na(up.marker.ch) && length(up.marker.ch)){
			# if(up.marker.ch == marker.ch){
				# up.marker.pos <- marker.bp[up.marker,2]
				# }else{
				# up.marker.pos <- 0
				# }
			# }else{
				# up.marker.pos <- 0
				# }
		# if(!is.na(down.marker.ch) && length(down.marker.ch) > 0){
			# if(down.marker.ch == marker.ch){
				# down.marker.pos <- marker.bp[down.marker,2]
				# }else{
				# down.marker.pos <- max.chr.len*1000000
				# }
			# }else{
				# down.marker.pos <- max.chr.len*1000000
			# }
			
		# result <- c(marker.ch, up.marker.pos, down.marker.pos)
		result <- c(marker.ch, min.pos, max.pos)
		names(result) <- c("chromosome", "min.boundary", "max.boundary")
		return(result)
		}
	
	#if we are only looking for specific blocks, trim down the list
	if(!is.null(block.num)){
		new.block.list <- vector(mode = "list", length = length(block.num))
		for(i in 1:length(block.num)){
			new.block.list[[i]] <- blocks[[block.num[i]]]
			}
		names(new.block.list) <- names(blocks)[block.num]
		blocks <- new.block.list
		}
	
	block.coords <- t(sapply(blocks, get.coord))
	write.table(block.coords, output.file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
	
	
	# min.table <- cbind(paste(rownames(block.coords), ".min", sep = ""), block.coords[,1], block.coords[,2])
	# max.table <- cbind(paste(rownames(block.coords), ".max", sep = ""), block.coords[,1], block.coords[,3])
	# final.table <- rbind(min.table, max.table)
	
	# write.table(final.table, output.file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	
}