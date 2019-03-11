#this script picks peak markers in 2 dimensions based
#on interaction data
#keep in mind that m12 is the effect of 2 on 1, not 
#the reverse

peak.markers.2d <- function(data.obj, chr = NULL, effect.drop = 1, pop.out.to.nearest.marker = TRUE, plot.peaks = TRUE, zlim = NULL){

library(igraph)
library(lattice)


if(!is.null(chr) && length(chr) < 2){
	stop("chr must have at least 2 elements.")
	}

if(!is.null(zlim) && length(zlim) != 2){
	stop("zlim must be a vector of length 2 specifying the min and max of the z axis.")
	}

marker.names <- data.obj$marker.names
all.markers <- 1:length(data.obj$marker.names)
imp.marker.locale <- grep("loc", data.obj$marker.names)

if(length(imp.marker.locale) > 0){
	true.marker.locale <- all.markers[-imp.marker.locale]
	}else{
	true.marker.locale <- all.markers
	}


pair.effects <- data.obj$var.to.var.influences
m12.std.effects <- abs(pair.effects[,3]/pair.effects[,4])
m21.std.effects <- abs(pair.effects[,5]/pair.effects[,6])
overall.min <- min(c(m12.std.effects, m21.std.effects))
overall.max <- max(c(m12.std.effects, m21.std.effects))

get.marker.chr <- function(marker.names){
	marker.chr <- data.obj$chromosome[match(marker.names, colnames(data.obj$geno))]
	return(marker.chr)
	}

all.markers.tested <- pair.effects[,1:2]
all.chr.tested <- apply(all.markers.tested, 2, get.marker.chr)




if(is.null(chr)){
	chr <- unique(as.vector(all.chr.tested))
	chr.pairs <- unique(all.chr.tested)
	}else{

	#check to see if the specified chromosomes are included in the results
	chr.tested <- chr %in% all.chr.tested
	if(length(which(chr.tested)) < length(chr)){
		cat("There were no markers tested on the following chromosomes:", chr[which(!chr.tested)], sep = "\n")
		return(invisible(NULL))
		}
			
	chr.locale1 <- which(all.chr.tested[,1] %in% chr); chr.locale2 <- which(all.chr.tested[,2] %in% chr)
	only.chr <- intersect(chr.locale1, chr.locale2)
	chr.pairs <- unique(all.chr.tested[only.chr,])
	}

	



	pop.to.markers <- function(marker.labels){

		final.true.markers <- rep(NA, 2)
		
		for(i in 1:length(marker.labels)){
			#find the chromosome of the marker
			marker.ch <- data.obj$chromosome[which(colnames(data.obj$geno) == marker.labels[i])]
			all.chr.markers <- colnames(data.obj$geno)[which(data.obj$chromosome == marker.ch)]
			true.chr.markers <- intersect(colnames(data.obj$geno)[true.marker.locale], all.chr.markers)

			#figure out of the marker is a true or imputed marker
			true.marker <- length(which(true.chr.markers == marker.labels[i]))
		
			if(true.marker == 1 || pop.out.to.nearest.marker == FALSE){
				final.true.markers[i] <- marker.labels[i]
				}else{
				#locate the marker among all markers
				imp.marker.locale <- which(all.chr.markers == marker.labels[i])
				true.chr.marker.locale <- which(all.chr.markers %in% true.chr.markers)
				if(i == 1){
					final.true.markers[i] <- all.chr.markers[max(true.chr.marker.locale[which(true.chr.marker.locale < imp.marker.locale)])]
					}else{
					final.true.markers[i] <- all.chr.markers[min(true.chr.marker.locale[which(true.chr.marker.locale > imp.marker.locale)])]	
					}
				}
			}
			
		marker.locale <- which(colnames(data.obj$geno) %in% final.true.markers)		
		marker.range <- min(marker.locale):max(marker.locale)
		marker.names <- data.obj$marker.names[marker.range]
		marker.chr <- data.obj$chromosome[marker.range]
		marker.pos <- data.obj$marker.location[marker.range]
		result <- matrix(c(marker.names, marker.chr, marker.pos), ncol = 3, byrow = FALSE)
		rownames(result) <- colnames(data.obj$geno)[marker.range]
		colnames(result) <- c("marker", "chromosome", "position")
		return(result)
		}


		get.peak <- function(m.mat){
			peak <- max(m.mat, na.rm = TRUE)
			min.m <- peak - effect.drop
			m.adj <- m.mat
			m.adj[which(m.mat < min.m)] <- 0
			
			# quartz();myImagePlot(m.mat)
			# test.mat <- m.mat; test.mat[which(test.mat < min.m)] <- 0
			# quartz();myImagePlot(test.mat)			
			
			mat.lat <- graph.lattice(dimvector = dim(m.mat))
			
			V(mat.lat)$vals <- unlist(apply(m.mat, 2, list))
			# V(mat.lat)$vals <- unlist(apply(m.mat, 1, list))
			
			dummy.mat <- matrix(1, nrow = dim(m.mat)[1], ncol = dim(m.mat)[2])
			all.ind <- which(dummy.mat == 1, arr.ind = TRUE)
			v.names <- apply(all.ind, 1, function(x) paste(rownames(m.mat)[x[1]], colnames(m.mat)[x[2]], collapse = ","))
			V(mat.lat)$name <- v.names
			V(mat.lat)$rows <- all.ind[,1]
			V(mat.lat)$cols <- all.ind[,2]
			V(mat.lat)$colors <- 1:vcount(mat.lat)

			# pdf("big.graph.pdf", width = 20, height = 20)
			# plot(mat.lat, vertex.label = V(mat.lat)$name, vertex.size = 1, vertex.color = heat.colors(vcount(mat.lat)))
			# dev.off()
			
			low.verts <- which(V(mat.lat)$vals < min.m)
			if(length(low.verts) > 0){
				mat.lat <- delete.vertices(mat.lat, low.verts)
				}
			na.verts <- which(is.na(V(mat.lat)$vals))
			if(length(na.verts) > 0){
				mat.lat <- delete.vertices(mat.lat, na.verts)
				}
			
			all.clusters <- igraph::clusters(mat.lat)
			# quartz();plot(mat.lat, vertex.color = ((all.clusters$membership)+1), vertex.label = V(mat.lat)$name)
		
			
			#take the cluster with the largest maximum
			max.val.locale <- which(V(mat.lat)$vals == max(V(mat.lat)$vals))[1]
			
			peak.cluster <- all.clusters$membership[max.val.locale]
			node.id <- which(all.clusters$membership == peak.cluster)
			
			node.names <- V(mat.lat)$name[node.id]
			
			
			clust.row.min <- min(V(mat.lat)$rows[node.id]); clust.row.max <- max(V(mat.lat)$rows[node.id])
			clust.col.min <- min(V(mat.lat)$cols[node.id]); clust.col.max <- max(V(mat.lat)$cols[node.id])
			
			markers.ch1 <- rownames(m.mat)[c(clust.row.min, clust.row.max)]
			markers.ch2 <- colnames(m.mat)[c(clust.col.min, clust.col.max)]
			
			
			ch1.range <- pop.to.markers(marker.labels = markers.ch1)
			ch2.range <- pop.to.markers(marker.labels = markers.ch2)
			
			ch1.range <- cbind(ch1.range, rep("Source", dim(ch1.range)[1]))
			ch2.range <- cbind(ch2.range, rep("Target", dim(ch2.range)[1]))
			
			results <- rbind(ch1.range, ch2.range)	
			
			return(results)			
			}
			

	plot.all <- function(m.mat, m.peaks, pair.num, m12.or.m21, debugging = FALSE){

		if(m12.or.m21 == "m12"){
			first.ch <- 2; second.ch = 1
			target.names <- list("Target", "Source")			
			}else{
			first.ch <- 1; second.ch = 2
			target.names <- list("Source", "Target")

			}
		
		z.minmax <- c(z.min, z.max)
		color.int <- (z.max - z.min)/25

		
		if(length(which(!is.na(m.mat))) < 2){
				for(i in 1:4){
				plot.new()
				plot.window(xlim = c(0,1), ylim = c(0,1))
				text(0.5, 0.5, paste("Not enough values to plot for Chr", chr.pairs[pair.num,first.ch], " -> ", "Chr", chr.pairs[pair.num,second.ch], sep = ""))
				}
			}else{
				#plot the 3D wire frame	
				if(m12.or.m21 == "m21"){
				print(wireframe(rotate.mat(rotate.mat(m.mat)), drape = TRUE, colorkey = list(at = seq(z.min, z.max, color.int)), xlab = paste("Chromosome", chr.pairs[pair.num,1]), ylab = paste("Chromosome", chr.pairs[pair.num,2]), zlab = paste("Ch", chr.pairs[pair.num,first.ch], " -> ", "Ch", chr.pairs[pair.num,second.ch], sep = ""), main = paste(m12.or.m21, ": Chr", chr.pairs[pair.num,first.ch], " -> ", "Chr", chr.pairs[pair.num,second.ch], sep = ""), zlim = z.minmax))
				}else{
				print(wireframe(rotate.mat(m.mat), drape = TRUE, colorkey = list(at = seq(z.min, z.max, color.int)), xlab = paste("Chromosome", chr.pairs[pair.num,1]), ylab = paste("Chromosome", chr.pairs[pair.num,2]), zlab = paste("Ch", chr.pairs[pair.num,second.ch], " -> ", "Ch", chr.pairs[pair.num,first.ch], sep = ""), main = paste(m12.or.m21, ": Chr", chr.pairs[pair.num,first.ch], " -> ", "Chr", chr.pairs[pair.num,second.ch], sep = ""), zlim = z.minmax))					
				}
				
				#plot the heatmap
				if(debugging){quartz()}
				myImagePlot(m.mat, main = paste(m12.or.m21, ": Chr", chr.pairs[pair.num,first.ch], " -> ", "Chr", chr.pairs[pair.num,second.ch], sep = ""), cex.main = 3)
		
				#plot the filtered heatmap
				if(debugging){quartz()}
				peak.mat <- m.mat; peak.mat[which(m.mat < (max(m.mat)-effect.drop))] <- 0
				myImagePlot(peak.mat, main = paste(m12.or.m21, " Peaks: Chr", chr.pairs[pair.num,second.ch], "->", "Chr", chr.pairs[pair.num,first.ch], sep = " "), cex.main = 3)

				#plot the actual region selected		
				selected.region  <- matrix(0, nrow = dim(m.mat)[1], ncol = dim(m.mat)[2])
				rownames(selected.region) <- rownames(m.mat); colnames(selected.region) <- colnames(m.mat)
				row.col <- lapply(target.names, function(x) which(m.peaks[,4] == x))
				selected.region[as.character(names(row.col[[first.ch]])),as.character(names(row.col[[second.ch]]))] <- 1
				if(debugging){quartz()}
				if(length(which(selected.region == 0)) == 0){
					selected.region <- jitter(selected.region, factor = 0.1)
					}	
					myImagePlot(selected.region, main = paste(m12.or.m21, " Peaks: Chr", chr.pairs[pair.num,second.ch], "->", "Chr", chr.pairs[pair.num,first.ch], sep = " "), min.x = 0, max.x = max(selected.region), cex.main = 3)
					# persp3d(m.mat, xlab = paste("Chromosome", chr.pairs[pair.num,first.ch]), ylab = paste("Chromosome", chr.pairs[pair.num,second.ch]), zlab = paste("Ch", chr.pairs[pair.num,first.ch], " -> ", "Ch", chr.pairs[pair.num,second.ch], sep = ""), col = heat.colors(100))
				}
			
		}

	all.interaction.ranges <- NULL

	
	#for each chromosome pair, find all the interaction effect sizes
	if(plot.peaks){
		pdf(paste("Peaks.2D.drop",effect.drop, ".pdf", sep = ""), width = 12, height = 12)
		layout(matrix(1:4, ncol = 2, byrow = TRUE))
		}
	for(i in 1:length(chr.pairs[,1])){
		# print(i)
		report.progress(i, dim(chr.pairs)[1])
		chr1.markers <- colnames(data.obj$geno)[which(data.obj$chromosome == chr.pairs[i,1])]
		chr2.markers <- colnames(data.obj$geno)[which(data.obj$chromosome == chr.pairs[i,2])]
		
		chr1.locale <- which(pair.effects[,1] %in% as.numeric(chr1.markers))
		chr2.locale <- which(pair.effects[,2] %in% as.numeric(chr2.markers))
		chr.int.locale <- intersect(chr1.locale, chr2.locale)
		
		#make the m12 and m21 matrices with source markers in rows and target markers in columns				
		m12.mat <- matrix(NA, nrow = length(chr2.markers), ncol = length(chr1.markers))
		m21.mat <- matrix(NA, nrow = length(chr1.markers), ncol = length(chr2.markers))


			colnames(m12.mat) <- rownames(m21.mat) <- chr1.markers
			colnames(m21.mat) <- rownames(m12.mat) <- chr2.markers

			for(j in 1:length(chr.int.locale)){
				m12.mat[as.character(pair.effects[chr.int.locale[j],]['marker2']), as.character(pair.effects[chr.int.locale[j],]["marker1"])] <- m12.std.effects[chr.int.locale[j]]
				m21.mat[as.character(pair.effects[chr.int.locale[j],]["marker1"]), as.character(pair.effects[chr.int.locale[j],]["marker2"])] <- m21.std.effects[chr.int.locale[j]]
				}
				
			# m.peaks <- get.peak(m.mat = m.mat)
			m12.peaks <- get.peak(m.mat = m12.mat)
			m21.peaks <- get.peak(m.mat = m21.mat)
	
			# m.peaks <- cbind(m.peaks, rep(paste(chr.pairs[i,2], "->", chr.pairs[i,1]), dim(m.peaks)[1]))			
			m12.peaks <- cbind(m12.peaks, rep(paste(chr.pairs[i,2], "->", chr.pairs[i,1]), dim(m12.peaks)[1]))
			m21.peaks <- cbind(m21.peaks, rep(paste(chr.pairs[i,1], "->", chr.pairs[i,2]), dim(m21.peaks)[1]))
						
			# all.interaction.ranges <- rbind(all.interaction.ranges, m.peaks)
			final.ranges <- rbind(m12.peaks, m21.peaks)
			all.interaction.ranges <- rbind(all.interaction.ranges, final.ranges)
			
				if(plot.peaks){	
					plot.all(m.mat = m12.mat, m.peaks = m12.peaks, pair.num = i, m12.or.m21 = "m12", debugging = FALSE)
					plot.all(m.mat = m21.mat, m.peaks = m21.peaks, pair.num = i, m12.or.m21 = "m21", debugging = FALSE)
					}
					

			}
		if(plot.peaks){
			dev.off()
			}
		
	colnames(all.interaction.ranges) <- c("marker", "chromosome", "position", "source.target", "direction")
	write.table(all.interaction.ranges, file = paste("Peaks.2D.drop",effect.drop, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	
}