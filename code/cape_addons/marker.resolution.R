#This script calculates the resolution of a cross
#using correlations between markers.
#it counts the number of LD blocks at the given
#r2.threshold and plots the sizes of the blocks
#along the genome

marker.resolution <- function(data.obj, r2.thresh = 0.8, pdf.label = "Genetic.Resolution.pdf"){
	
	old.r2.thresh <- data.obj$r2.thresh
	if(is.null(old.r2.thresh) || old.r2.thresh != r2.thresh){
		cat("running linkage.blocks() to determine block structure...\n")
		data.obj <- all.linkage.blocks(data.obj, r2.thresh = r2.thresh)
		}

	link.blocks <- data.obj$linkage.blocks

	
	layout.mat <- matrix(c(1,1,1,2,3,4,5,5,5,6,6,6), ncol = 3, byrow = TRUE)
	
	
	pdf(pdf.label, width = 12, height = 12)
	layout(layout.mat, heights = c(1,1,0.5,1))
	
	marker.cor <- cor(data.obj$geno, use = "complete")
	diag(marker.cor) <- NA
	#==============================================================
	# plot 1
	# boxplot of correlations for each marker with all other 
	# markers
	#==============================================================
	boxplot(abs(marker.cor), ylim = c(0, 1), axes = FALSE, main = "Correlations By Marker")
	abline(h = r2.thresh)
	axis(2)
	axis(1, at = 1:dim(marker.cor)[1], labels = FALSE)
	par(xpd = TRUE)
	text(x = 1:dim(marker.cor)[1], y = rep(-0.09, dim(marker.cor)[1]), data.obj$marker.names, srt = 90, cex = 0.6, adj = 1)
	par(xpd = FALSE)
		
	#==============================================================
	# plot 2
	# histogram of all pairwise correlations between markers
	#==============================================================	
	marker.cor[lower.tri(marker.cor, diag = TRUE)] <- NA
	
	hist(abs(as.vector(marker.cor)), main = "All Pairwise Correlations", freq = FALSE, xlab = "Correlation", xlim = c(0, 1))
	abline(v = r2.thresh)

	#find the average distance between all markers on the same chromosome
	chr <- unique(data.obj$chromosome)
	num.chr <- length(chr)
	distances <- rep(NA, (dim(marker.cor)[1]-1))
	for(i in 1:(dim(marker.cor)[1]-1)){
		if(data.obj$chromosome[i] == data.obj$chromosome[(i+1)]){
			distances[i] <- abs(data.obj$marker.location[i] - data.obj$marker.location[(i+1)])
			}
		}
		
	#==============================================================
	# plot 3
	# histogram of distances between adjacent markers
	#==============================================================		
	hist(distances, main = "Distances Between Adjacent Markers", xlab = "Distance")
			
	num.blocks <- length(link.blocks)
	markers.per.block <- sapply(link.blocks, length)
	large.blocks <- which(markers.per.block > 1)



	#for plotting, calculate beginning and end x 
	#coordinates of each chromosome
	chr.x <- matrix(NA, ncol = 2, nrow = num.chr)
	rownames(chr.x) <- chr
	chr.x[,1] <- (0:(num.chr-1))+(num.chr*0.005)
	chr.x[,2] <- (1:num.chr)-(num.chr*0.005)

	#calculate the size of the blocks 
	block.position <- matrix(NA, nrow = num.blocks, ncol = 2)
	all.chr <- rep(NA, num.blocks)
	colnames(block.position) <- c("start", "stop")
	for(i in 1:num.blocks){
		start.block <- link.blocks[[i]][1]
		end.block <- link.blocks[[i]][length(link.blocks[[i]])]
		
		marker1.locale <- which(colnames(data.obj$geno) == start.block)
		marker2.locale <- which(colnames(data.obj$geno) == end.block)
		
		block.position[i,1] <- data.obj$marker.location[marker1.locale]
		block.position[i,2] <- data.obj$marker.location[marker2.locale]

		all.chr[i] <- data.obj$chromosome[marker1.locale]
		}
	
	relative.block.position <- matrix(NA, nrow = num.blocks, ncol = 2)
	for(i in 1:length(chr)){
		chr.locale <- which(all.chr == chr[i])
		if(length(chr.locale)){
			chr.max <- max(block.position[chr.locale,])
			relative.block.position[chr.locale,] <- ((block.position[chr.locale,]/chr.max)*(chr.x[1,2]-chr.x[1,1]))+chr.x[i,1]
			bad.vals <- which(!is.finite(relative.block.position[chr.locale,]))
			if(length(bad.vals > 0)){
				relative.block.position[chr.locale,bad.vals] <- chr.x[i,1]
				}
			}
		}
			

	#==============================================================
	# plot 4
	# histogram of sizes of blocks with more than one marker
	#==============================================================			
	block.sizes <- abs(block.position[large.blocks,2] - block.position[large.blocks,1])	
	hist(block.sizes, xlab = "Linkage Block Size", main = paste("Size of Blocks With >1 Marker\nR2 = ", r2.thresh))
	

	#==============================================================
	# plot 5
	# text describing numerical results
	#==============================================================		
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	text(0.2, 0.5, paste(length(link.blocks), "blocks for", dim(marker.cor)[1], "markers at R2 = ", r2.thresh), cex = 2)
	text(0.7, 0.5, paste(length(large.blocks), "blocks with more than 1 marker."), cex = 2)


	# dev.new(width = 12, height = 3)
	plot.new()
	plot.window(xlim = c(0,length(chr)), ylim = c(0,1))
	par(mar = c(0,0,0,0))
	for(i in 1:length(chr.x[,1])){
		segments(chr.x[i,1], 0.5, chr.x[i,2], 0.5, lwd = 5)
		text(x = mean(chr.x[i,]), y = 0.3, chr[i], cex = 1.5)
		}
	for(i in 1:dim(relative.block.position)[1]){
		polygon(x = c(relative.block.position[i,], rev(relative.block.position[i,])), y = c(0.5,0.5,0.7,0.7), col = "red")
		}
	mtext("Block Relative Positions and Sizes", lin = -3)
	
	dev.off()
	
	
	
}