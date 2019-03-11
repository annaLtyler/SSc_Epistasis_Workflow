#This function takes the results stored in data.obj$p.adjusted
#and plots a network.
#In this version, the main effects are embeded in the nodes

plot.network.pheno.pie <- function(data.obj, p.or.q = 0.05,  collapsed.net = FALSE, layout.matrix = NULL, node.radius = 1, label.nodes = TRUE, label.offset = 0, label.cex = 1, legend.radius = 1, legend.cex = 1, arrow.offset = 1, edge.lwd = 2, singleton.position = "bottom", layout.type = c("auto", "kamada.kawai", "circle")){
	
	
	require(igraph)
		
	iArrows <- igraph:::igraph.Arrows		

	pheno.tables <- data.obj$max.var.to.pheno.influence
	phenotypes <- names(pheno.tables)	


	plot.net.edges <- function(net, net.layout, lwd = 1, edge.col = "gray", arrow.offset = 1){
	
		edge.list <- get.edgelist(net, names = FALSE)
	
		if(length(col) == 1){
			edge.col = rep(edge.col, dim(edge.list)[1])
			}


		for(i in 1:length(edge.list[,1])){
			xy.start = c(net.layout[edge.list[i,1],1], net.layout[edge.list[i,1],2])
			xy.end <- c(net.layout[edge.list[i,2],1], net.layout[edge.list[i,2],2])
			
			if(edge.list[i,1] != edge.list[i,2]){
				alpha <- 0
				final.xy <- alpha*xy.start + (1-alpha)*(xy.end)
				dist.to.center <- dist(matrix(c(xy.end, final.xy), nrow = 2, byrow = TRUE), method = "euclidean")
				while(dist.to.center < arrow.offset){
					alpha <- alpha + 0.01
					final.xy <- alpha*xy.start + (1-alpha)*(xy.end)
					dist.to.center <- dist(matrix(c(xy.end, final.xy), nrow = 2, byrow = TRUE), method = "euclidean")	
					}
				# cat(round(xy.start, 2), " : ", round(final.xy, 2), "\n")
				arrows(x0 = xy.start[1], x1 = final.xy[1], y0 = xy.start[2], y1 = final.xy[2], lwd = lwd, col = edge.col[i], length = 0.2)
				}else{
				xy.start <- xy.start+node.radius*0.9; xy.end <- xy.end-node.radius*0.9
				iArrows(x1 = xy.start[1], y1 = xy.start[2], x2 = xy.end[1], y2 = xy.end[2], h.lwd=lwd, sh.lwd=lwd, sh.col=edge.col, curve = 2, width=1, size=0.7)	
				}
			}
		
		}

	#This function plots nodes of a graph as polygons where each polygon
	#is subdivided into polygons each with a different color
	#This is for adding information to nodes about which phenotypes they
	#have main effects on
	plot.net.point <- function(x, y, cols = c("green", "green", "red"), edge.len = 0.025, node.label = NULL, label.offset = 0, label.cex = label.cex){
		
		draw.pie(x = x, y = y, radius = node.radius, cols = cols, add = TRUE)
		
		
		# min.x <- (x-edge.len); max.x <- (x+edge.len)
		# min.y <- (y-edge.len); max.y <- (y+edge.len)
		
		# even.x <- segment.region(region.min = min.x, region.max = max.x, num.points = (length(cols)+1), alignment = "ends")
		
		# for(i in 1:length(cols)){
			# polygon(x = c(even.x[i], even.x[(i+1)], even.x[(i+1)], even.x[i]), y = c(min.y, min.y, max.y, max.y), col = cols[i])
			# }
		
		if(!is.null(node.label)){
			text.x <- x+label.offset; text.y = y+label.offset
			text(x = text.x, y = text.y, labels = node.label, cex = label.cex)
			}
		
		}





	if(collapsed.net){
		net.data <- data.obj$collapsed.net
		if(is.null(net.data)){
			stop("get.network needs to be run with collapsed = TRUE before plotting the collapsed network.")
			}
		}else{
		net.data <- data.obj$full.net
		if(is.null(net.data)){
			stop("get.network needs to be run with collapsed = FALSE before plotting the collapsed network.")
			}			
		}
		
		#convert this into an edge list to be compatible with the uncollapsed network
		sig.locale <- which(net.data != 0, arr.ind = TRUE)
			
			if(length(sig.locale) == 0){
				plot.new()
				plot.window(xlim = c(0,1), ylim = c(0,1))
				text(0.5, 0.5, "No Significant Interactions")
				invisible()
				}
		edge.vals <- net.data[which(net.data != 0)]
		sig.edges <- cbind(sig.locale, edge.vals)
		colnames(sig.edges) <- c("Source", "Target", "Effect")
		block.names <- matrix(NA, ncol = 2, nrow = dim(sig.edges)[1])
		block.names[,1] <- rownames(net.data)[sig.edges[,1]]
		block.names[,2] <- colnames(net.data)[sig.edges[,2]]
		sig.edges[,1:2] <- block.names
					
	
		get.node.cols <- function(node.name, phenotypes, main.effects){
			node.locale <- which(main.effects[,"Source"] == node.name)
			node.cols <- rep("gray", length(phenotypes))
			
			node.effects <- main.effects[node.locale,c("Target", "Effect"),drop=FALSE]
			
			if(length(node.effects) > 0){
				for(i in 1:length(phenotypes)){
					pheno.locale <- which(node.effects[,"Target"] == phenotypes[i])
					if(length(pheno.locale) > 0){
						if(as.numeric(node.effects[pheno.locale,"Effect"]) > 0){node.cols[i] <- "green"}
						if(as.numeric(node.effects[pheno.locale,"Effect"]) < 0){node.cols[i] <- "red"}
						}
					}
				}
			return(node.cols)
			}
			
						
		main.effects.locale <- which(sig.edges[,2] %in% phenotypes)
		main.effects <- sig.edges[main.effects.locale,]
		interactions <- sig.edges[-main.effects.locale,]
		
		all.interacting.nodes <- unique(as.vector(c(interactions[,1], interactions[,2])))
		all.main.effect.nodes <- unique(as.vector(c(main.effects[,1], main.effects[,2])))
		non.interacting <- setdiff(all.main.effect.nodes, all.interacting.nodes)
		unordered.singletons <- setdiff(non.interacting, phenotypes)
		singletons <- rownames(net.data)[sort(match(unordered.singletons, rownames(net.data)))]
		
		edgelist <- matrix(c(as.vector(interactions[,1]), as.vector(interactions[,2])), ncol = 2, byrow = FALSE)
	
		net <- graph.edgelist(edgelist)
		vertex.names <- V(net)$name
		# if(collapsed.net){
			# vertex.names <- V(net)$name
			# }else{
			# split.names <- apply(matrix(V(net)$name, ncol = 1), 1, function(x) strsplit(x, "_"))
			# marker.loc <- sapply(split.names, function(x) x[[1]][2])
			# all.split <- lapply(split.names, function(x) strsplit(x[[1]][1], "Chr"))
			# marker.chr <- sapply(all.split, function(x) x[[1]][2])
			# vertex.names <- rep(NA, length(marker.chr))
			# for(j in 1:length(marker.chr)){
				# vertex.names[j] <- data.obj$marker.names[which(data.obj$chromosome == marker.chr[j])[marker.loc[j]]]
				# }
			# na.locale <- which(is.na(vertex.names))
			# vertex.names[na.locale] <- V(net)$name[na.locale]
			# }
			
		if(!label.nodes){vertex.names <- NULL}
		
		edge.col <- rep(NA, ecount(net))
		edge.col[which(as.numeric(interactions[,"Effect"]) < 0)] <- "red"
		edge.col[which(as.numeric(interactions[,"Effect"]) > 0)] <- "green"
				
		if(class(layout.matrix) != "matrix" && !is.null(layout.matrix)){
			if(layout.matrix == "manual"){
				tkp.id <- tkplot(net)
				done <- readline(prompt = "Press return when ready:\n")
				layout.matrix <- tkplot.getcoords(tkp.id)
				tkplot.close(tkp.id)
				}
			}


		if(is.null(layout.matrix)){
			if(length(layout.type) > 1 || layout.type == "auto"){
				layout.type = "auto"
				coord.matrix <- layout.auto(net)
				}
			if(layout.type == "kamada.kawai"){
				coord.matrix <- layout.kamada.kawai(net)	
				}
			if(layout.type == "circle"){
				coord.matrix <- layout.circle(net)
				# coord.matrix <- layout.svd(net)
				# coord.matrix <- layout.lgl(net)
				# coord.matrix <- layout.graphopt(net)
				# coord.matrix <- layout.fruchterman.reingold.grid(net)
				}
			}else{
			coord.matrix <- layout.matrix	
			}

			#scale the coordinates so that 
			#the axes go from 0 to plot.dim
			plot.dim <- 150
			if(min(coord.matrix) < 0){
				abs.coord <- coord.matrix + abs(min(coord.matrix))
				}else{
				abs.coord <- coord.matrix	
				}
				
			max.coord <- max(abs.coord)
			scaling.factor <- plot.dim/max.coord
			scaled.coord <- abs.coord*scaling.factor
			
			#center the network in the plotting area
			xy.offset <- (plot.dim/2) - apply(scaled.coord, 2, mean)
			xy.min <- apply(scaled.coord, 2, min)
			xy.max <- apply(scaled.coord, 2, max)
			#the x and y coordinates can be shifted in the direction
			#of the xy.offset, but not more than the extremes of the
			#current coordinates
			xy.shift <- rep(NA, 2)
			for(i in 1:2){
				if(xy.offset[i] < 0){
					xy.shift[i] <- min(c(abs(xy.offset[i]), xy.min[i]))*-1
					}else{
					xy.shift[i] <- min(c(abs(xy.offset[i]), xy.max[i]))	
					}
				}
						
			centered.coord <- matrix(NA, nrow = dim(scaled.coord)[1], ncol = 2)
			for(i in 1:2){
				centered.coord[,i] <- scaled.coord[,i]+xy.shift[i]
				}
			
			#if there are singletons, make a coordinate matrix for these too.
			#place them along the bottom edge of the network
			if(length(singletons) > 0){
				singleton.coords <- matrix(NA, ncol = 2, nrow = length(singletons))
				if(singleton.position == "bottom"){
					singleton.coords[,1] <- segment.region(min(centered.coord[,1]), max(centered.coord[,1]), length(singletons), "ends")
					singleton.coords[,2] <- rep((0-(plot.dim*0.03)), length(singletons))
					}
				if(singleton.position == "top"){
					singleton.coords[,1] <- segment.region(min(centered.coord[,1]), max(centered.coord[,1]), length(singletons), "ends")
					singleton.coords[,2] <- rep(max(centered.coord[,2])*1.07, length(singletons))					
					}
				if(singleton.position == "left"){
					singleton.coords[,1] <- rep(-1, length(singletons))
					singleton.coords[,2] <- segment.region(min(centered.coord[,1]), max(centered.coord[,1]), length(singletons), "center")
					centered.coord[,1] <- centered.coord[,1] + plot.dim*0.05
					}
				if(singleton.position == "right"){
					singleton.coords[,1] <- rep(max(centered.coord[,1])*1.07, length(singletons))
					singleton.coords[,2] <- segment.region(min(centered.coord[,1]), max(centered.coord[,1]), length(singletons), "ends")
					centered.coord[,1] <- (centered.coord[,1] - plot.dim*0.05)
					}
				}
			
			# quartz(width = 8, height = 8)
			# par(mfrow = c(2,2))
			# plot(coord.matrix[,1], coord.matrix[,2], pch = 16, col = "blue")
			# plot(abs.coord[,1], abs.coord[,2], pch = 16, col = "blue")
			# plot(scaled.coord[,1], scaled.coord[,2], pch = 16, col = "blue", xlim = c(0,plot.dim), ylim = c(0,plot.dim))
			# plot(centered.coord[,1], centered.coord[,2], pch = 16, col = "blue", xlim = c(0,plot.dim), ylim = c(0,plot.dim))

		
		if(!is.null(layout.matrix) & length(layout.matrix[,1]) != vcount(net)){
			stop("The layout matrix does not have the same number of nodes as the network.")
			}


		border = plot.dim*0.01
		if(label.nodes){
			singleton.labels <- singletons
			}else{
			singleton.labels <- rep(NA, length(singletons))
			}
		
	
		plot.new()
		plot.window(xlim = c((0-border), (plot.dim+border)), ylim = c((0-border), (plot.dim+border)))
		plot.net.edges(net = net, net.layout = centered.coord, lwd = edge.lwd, edge.col = edge.col, arrow.offset = arrow.offset)
		for(i in 1:length(V(net)$name)){
			plot.net.point(x = centered.coord[i,1], y = centered.coord[i,2], cols = get.node.cols(node.name = V(net)$name[i], phenotypes = phenotypes, main.effects = main.effects), edge.len = node.box.length, node.label = vertex.names[i], label.offset = label.offset, label.cex = label.cex)
			}
		if(length(singletons) > 0){
			for(i in 1:length(singletons)){
				plot.net.point(x = singleton.coords[i,1], y = singleton.coords[i,2], cols = get.node.cols(node.name = singletons[i], phenotypes = phenotypes, main.effects = main.effects), edge.len = node.box.length, node.label = singleton.labels[i], label.offset = label.offset, label.cex = label.cex)		
				}
			}
		draw.pie(x = 0, y = plot.dim, radius = legend.radius, cols = rep("gray", length(phenotypes)), labels = phenotypes, label.cex = legend.cex)
		
	
		results <- list(net, centered.coord)
		names(results) <- c("net", "layout.matrix")
		invisible(results)

	
	}
