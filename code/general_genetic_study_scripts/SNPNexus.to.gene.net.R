#This snp.gene.table used here is output from SNPNexus
# snp.gene.table <- as.matrix(read.table("~/Documents/Data/Scleroderma/Results/Dominant_lung_pheno_2ET/SSc_snpnexus_11708/ncsnp_11708.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
# require(igraph)
# require(RColorBrewer)


SNPNexus.to.gene.net <- function(data.obj, snp.gene.table, p.or.q = 0.05, vertex.size = 50, label.x.shift = 0, label.y.shift = 0, cex.label = 1){
	
	# display.brewer.all()
	var.inf <- writeVariantInfluences(data.obj, p.or.q, include.main.effects = FALSE, write.file = FALSE)
	edge.list <- var.inf[,c(1,4)]	
	edge.list[,1] <- sapply(strsplit(edge.list[,1], "_"), function(x) x[1])
	edge.list[,2] <- sapply(strsplit(edge.list[,2], "_"), function(x) x[1])	
	
	#======================================================
	# internal functions
	#======================================================
	
	snp.gene <- function(snp){
		snp.locale <- which(snp.gene.table[,2] == snp)
		if(length(snp.locale) == 0){
			return(rep(snp, 4))
			}
		
		overlap.gene <- snp.gene.table[snp.locale,"Overlapped.Gene"]
		if(overlap.gene != ""){
			gene.info <- c(overlap.gene, "overlapping", NA, snp.gene.table[snp.locale,"Type"])
			}else{
			upstream.dist <- as.numeric(snp.gene.table[snp.locale,"Distance.to.Nearest.Upstream.Gene"])
			downstream.dist <- as.numeric(snp.gene.table[snp.locale,"Distance.to.Nearest.Downstream.Gene"])
			
			if(upstream.dist < downstream.dist){
				gene.info <- c(snp.gene.table[snp.locale,"Nearest.Upstream.Gene"], "upstream", upstream.dist, snp.gene.table[snp.locale,"Type.of.Nearest.Upstream.Gene"])
				}else{
				gene.info <- c(snp.gene.table[snp.locale,"Nearest.Downstream.Gene"], "downstream", upstream.dist, snp.gene.table[snp.locale,"Type.of.Nearest.Downstream.Gene"])	
				}
			}
		return(gene.info)
		}
		
		
	#from shapes example
	mystar <- function(coords, v=NULL, params) {
	  vertex.color <- params("vertex", "color")
	  if (length(vertex.color) != 1 && !is.null(v)) {
	    vertex.color <- vertex.color[v]
	  }
	  vertex.size  <- 1/200 * params("vertex", "size")
	  if (length(vertex.size) != 1 && !is.null(v)) {
	    vertex.size <- vertex.size[v]
	  }
	  norays <- params("vertex", "norays")
	  if (length(norays) != 1 && !is.null(v)) {
	    norays <- norays[v]
	  }
	
	  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
	         FUN=function(x, y, bg, size, nor) {
	           symbols(x=x, y=y, bg=bg,
	                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
	                   add=TRUE, inches=FALSE)
	         })
		}


	assign.colors <- function(values, color.table){
		cols <- rep("gray", length(values))
		for(i in 1:nrow(color.table)){
			cols[which(values == color.table[i,1])] <- color.table[i,2]
			}
		return(cols)
		}
	#======================================================
	
	snp.net <- graph_from_edgelist(edge.list)
	edge.weight <- as.numeric(var.inf[,"Effect"])
	pos.col <- get.color("brown")[3]
	neg.col <- get.color("blue")[3]
	edge.color <- rep(pos.col, ecount(snp.net))
	edge.color[which(edge.weight < 0)] <- neg.col
	
	net.layout <- layout_with_kk(snp.net)

	# plot(snp.net)
	chr.cols <- c(brewer.pal(8, "Set2"), brewer.pal(9, "Set1"))
	u_chr <- unique(c(var.inf[,2], var.inf[,5]))
	chr.table <- cbind(u_chr, chr.cols[1:length(u_chr)])
	chr.table <- chr.table[order(as.numeric(chr.table[,1])),]
	
	marker.chr <- sapply(V(snp.net)$name, function(x) get.marker.chr(data.obj, x))
	chr.col <- assign.colors(marker.chr, chr.table)
	# plot(snp.net, vertex.color = chr.col, layout = net.layout)
	# legend("topleft", legend = chr.table[,1], fill = chr.table[,2], title = "Chromosome")
	
	all.gene.info <- t(sapply(V(snp.net)$name, snp.gene))


	overlap.colors <- cbind(c("overlapping", "upstream", "downstream"), brewer.pal(3, "Set2"))
	cols <- assign.colors(values = all.gene.info[,2], overlap.colors)
	
		
	u_types <- unique(all.gene.info[,"Type"])
	num.stars <- length(u_types) - 2
	all.shapes <- c("circle", "square", rep("star", num.stars))
	sides <- c(1, 1, 3:(num.stars+2))
	type.shapes <- cbind(u_types, all.shapes, sides)
	
	# no clipping, edges will be below the vertices anyway
	add_shape("star", clip=shape_noclip, plot=mystar, parameters=list(vertex.norays=5))
	
	v.shape <- rep(NA, vcount(snp.net))
	v.sides <- rep(NA, vcount(snp.net))
	for(i in 1:nrow(type.shapes)){
		v.shape[which(all.gene.info[,"Type"] == type.shapes[i,1])] <- type.shapes[i,2]
		v.sides[which(all.gene.info[,"Type"] == type.shapes[i,1])] <- as.numeric(type.shapes[i,3])
		}

	par(xpd = TRUE, mar = c(2,2,2,4))	
	layout(matrix(c(2,3,4,1,1,1), ncol = 2), widths = c(0.4,1))
	xlim <- c(min(net.layout[,1]), max(net.layout[,1])*1.2)
	ylim <- c(min(net.layout[,2]), max(net.layout[,2]))
	plot.new()
	plot.window(xlim = xlim, ylim = ylim)
	plot(snp.net, vertex.size = vertex.size, vertex.label = NA, vertex.color = cols, vertex.shape = v.shape, vertex.norays = v.sides, layout = net.layout, edge.color = edge.color, edge.width = 3, rescale = FALSE, add = TRUE)
	
	v.labels <- apply(cbind(rownames(all.gene.info), all.gene.info[,1]), 1, function(x) paste(unique(x), collapse = "\n"))
	par(xpd = TRUE)
	text(net.layout[,1]+label.x.shift, net.layout[,2]+label.y.shift, labels = v.labels, adj = 0, cex = cex.label)
	par(xpd = FALSE)
	
	
	#plot the legend for the colors
	plot.new();plot.window(xlim = c(0,1), ylim = c(0,1))
	par(mar = c(0,0,0,0))
	legend(0, 0.5, fill = overlap.colors[,2], legend = overlap.colors[,1], cex = 1.5)
	
	#make a network to plot the legend shapes
	adj.mat <- diag(0, nrow = nrow(type.shapes))
	legend.net <- graph_from_adjacency_matrix(adj.mat)	
	V(legend.net)$name <- type.shapes[,1]
	

	x <- rep(-1, nrow(type.shapes))
	y <- segment.region(0, 0.9, nrow(type.shapes), alignment = "ends")
	layout.mat <- cbind(x,y)
	
	#plot the legend for the vertex shapes
	par(mar = c(0,0,0,0))
	plot(legend.net, layout = layout.mat, vertex.shape = type.shapes[,2], vertex.color = "white", vertex.norays = as.numeric(type.shapes[,3]), rescale = FALSE, vertex.label = NA, vertex.size = vertex.size/2)
	text(layout.mat[,1]+label.x.shift/2, layout.mat[,2], V(legend.net)$name, adj = 0, cex = 1.5)
	par(xpd = FALSE)
	
	#plot the legend for edge colors
	plot.new();plot.window(xlim = c(0,1), ylim = c(0,1))
	legend(0.1, 0.5, lty = 1, col = c(pos.col, neg.col), legend = c("Enhancing", "Suppressing"), cex = 1.5, lwd = 3)
	
	V(snp.net)$gene.name <- all.gene.info[,1]
	invisible(snp.net)
	
	}