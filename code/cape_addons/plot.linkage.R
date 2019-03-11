#This function plots correlation between markers
#on a chromosome by chromosome basis
#This script also marks regions of high linkage
#depth controls how many pairs are checked for correlation
# The script only performs correlation if the markers are
#within the distance specified. This can be either in base
#pairs or in the centimorgans depending on the 
#if depth is set to -1, the linkage is shown between all
#markers regardless of chromosome.

plot.linkage <- function(data.obj, mark.above.r = 0.7, depth = NULL, by.chromosome = FALSE, show.marker.labels = FALSE, show.chr.boundaries = TRUE, label.chr = TRUE, filename = "Linkage.pdf"){

	
	all.geno <- data.obj$geno
	ind.geno <- data.obj$geno.for.pairscan
	chr <- data.obj$chromosome
	pos <- data.obj$marker.location
	
	u_chr <- unique(chr)
	
	unique.markers <- colnames(data.obj$geno)
		
	unique.marker.locale <- which(colnames(data.obj$geno) %in% unique.markers)
	#get coordinates of the chromosome boundaries
	if(show.chr.boundaries){
		chromosomes <- data.obj$chromosome[sort(unique.marker.locale)]
		u_chr <- unique(chromosomes[which(!is.na(chromosomes))])
		chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
		chr.boundaries <- c(0, chr.boundaries)
		if(label.chr){
			chr.names <- unique(chromosomes)
			}else{
			chr.names <- NULL
			}
		}else{
		chr.boundaries <- NULL
		chr.names <- NULL
		}


	if(by.chromosome){
	pdf(filename)
	for(ch in 1:length(u_chr)){
		report.progress(ch, length(u_chr), 10)
		
		markers <- which(chr == u_chr[ch])
		num.markers <- length(markers)
		all.pairs <- pair.matrix(markers)
		num.pairs <- dim(all.pairs)[1]
		
		ind.markers <- match(colnames(ind.geno), colnames(all.geno)[markers])
		# to.mark <- ind.markers[which(!is.na(ind.markers))]
		
		if(num.pairs > 1){
			all.distances <- apply(all.pairs, 1, function(x) abs(pos[x[1]] - pos[x[2]]))
			if(is.null(depth)){
				depth <- max(all.distances)
				}
	
			within.limit <- which(all.distances <= depth)	
			cor.mat <- matrix(NA, num.markers, num.markers)
			colnames(cor.mat) <- rownames(cor.mat) <- markers
			for(i in 1:num.pairs){
				test <- try(cor(all.geno[,all.pairs[i,1]], all.geno[,all.pairs[i,2]], use = "complete.obs"), silent = TRUE)
				if(class(test) != "try-error"){
					cor.mat[as.character(all.pairs[i,1]), as.character(all.pairs[i,2])] <- cor(all.geno[,all.pairs[i,1]], all.geno[,all.pairs[i,2]], use = "complete.obs")
					}
				# cor.mat[as.character(all.pairs[i,2]), as.character(all.pairs[i,1])] <- cor(all.geno[,all.pairs[i,1]], all.geno[,all.pairs[i,2]], use = "complete.obs")
				}
		
			colnames(cor.mat) <- rownames(cor.mat) <- colnames(all.geno)[markers]
			# if(length(to.mark) > 0){
				# for(i in 1:length(to.mark)){
					# cor.mat[to.mark[i], to.mark[i]] <- 100 #use an arbitrary number that will otherwise never show up in the matrix
					# }				
				# }
			selected.coords <- which(rotate.mat(cor.mat) == 100, arr.ind = TRUE)
			# cor.mat[which(cor.mat == 100)] <- NA
			# myImagePlot(cor.mat, main = paste("Chromosome", u_chr[ch]), col = topo.colors(50), mark.coords = selected.coords, mark.col = "red", min.x = 0, max.x = 1)
			myImagePlot(cor.mat, main = paste("Chromosome", u_chr[ch]), mark.coords = selected.coords, mark.col = "red", min.x = 0, max.x = 1, col.split.point = mark.above.r, pos.col = "brown", neg.col = "blue")
			# legend("bottomright", legend = "Marker Linearly Independent of all Other Markers", col = "red", pch = 16)
			} #end case for if we have more than one pair to look at		
		
		} #end looping through chromosomes
	dev.off()
	}else{

		all.pairs <- pair.matrix(1:dim(all.geno)[2])
		num.pairs <- dim(all.pairs)[1]
	
		all.distances <- apply(all.pairs, 1, function(x) abs(pos[x[1]] - pos[x[2]]))
	
		if(!is.null(depth) && depth == -1){
			same.chrom <- rep(TRUE, length(all.pairs[,1]))
			depth <- max(all.distances)
			}else{
			same.chrom <- apply(all.pairs, 1, function(x) (chr[x[1]] == chr[x[2]]))
			if(is.null(depth)){
				depth <- max(all.distances)
				}
			}
		

		within.limit <- intersect(which(all.distances <= depth), which(same.chrom))
		cor.mat <- matrix(NA, dim(all.geno)[2], dim(all.geno)[2])

		for(i in 1:length(within.limit)){
			test <- try(cor(all.geno[,all.pairs[within.limit[i],1]], all.geno[,all.pairs[within.limit[i],2]], use = "complete.obs"), silent = TRUE)
			if(class(test) != "try-error"){
				cor.mat[all.pairs[within.limit[i],1], all.pairs[within.limit[i],2]] <- suppressWarnings(cor(all.geno[,all.pairs[within.limit[i],1]], all.geno[,all.pairs[within.limit[i],2]], use = "complete.obs"))
				}
			}
	
				
		colnames(cor.mat) <- rownames(cor.mat) <- data.obj$marker.names
		
		ind.markers <- match(colnames(ind.geno), colnames(all.geno))
		myImagePlot(cor.mat, main = "Correlations Between Markers on All Chromosomes", min.x = 0, max.x = 1, show.labels = show.marker.labels, chromosome.coordinates = chr.boundaries, chr.names = chr.names, col.split.point = mark.above.r, pos.col = "brown", neg.col = "blue")
		}
	

	
	invisible(cor.mat)


}
