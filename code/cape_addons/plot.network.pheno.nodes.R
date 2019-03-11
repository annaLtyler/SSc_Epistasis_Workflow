#This function takes the results stored in data.obj$p.adjusted
#and plots a network.

plot.network.pheno.nodes <- function(data.obj, p.or.q = 0.05, weight.edges = FALSE, separate.pheno = FALSE, highlight.node.edges = NULL, layout.matrix = NULL, collapsed.net = FALSE){
	
	
	
	require(igraph)
		
	pheno.tables <- data.obj$max.var.to.pheno.influence
	phenotypes <- names(pheno.tables)	

	if(collapsed.net){
		net.data <- data.obj$collapsed.net
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
				
		}else{
			
			net.data <- data.obj$var.to.var.p.val
							
			var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(net.data))))
			sig.locale <- which(net.data[,var.sig.col] <= p.or.q)

			if(length(sig.locale) == 0){
				plot.new()
				plot.window(xlim = c(0,1), ylim = c(0,1))
				text(0.5, 0.5, "No Significant Interactions")
				invisible()
				}

			#get the edges for variant to variant interactions
			sig.edges <- matrix(net.data[sig.locale,1:2], nrow = length(sig.locale))
			effect.size <- (net.data[sig.locale,"Effect"]/net.data[sig.locale,"SE"])
			sig.edges <- cbind(sig.edges, effect.size)
			colnames(sig.edges) <- c("Source", "Target", "Effect")

				#add any edges from the phenotype direct influences
				for(ph in 1:length(phenotypes)){
					pheno.stats <- pheno.tables[[ph]][,which(colnames(pheno.tables[[ph]]) != "t.stat")]
					p.adj.locale <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.stats))))
					sig.ph.edges <- which(as.numeric(pheno.stats[,p.adj.locale]) <= p.or.q)
					if(length(sig.ph.edges) > 0){
						for(r in sig.ph.edges){
							# table.row <- c(pheno.stats[r,1], names(pheno.tables)[ph], 
							# as.numeric(pheno.stats[r,grep("coef", colnames(pheno.stats))]),
							# as.numeric(pheno.stats[r,grep("se", colnames(pheno.stats))]),
							# as.numeric(pheno.stats[r,grep("\\|t.stat\\|", colnames(pheno.stats))]),
							# as.numeric(pheno.stats[r,grep("emp.p", colnames(pheno.stats))]),
							# as.numeric(pheno.stats[r,p.adj.locale]))							
							table.row <- c(pheno.stats[r,1], names(pheno.tables)[ph], 
							as.numeric(pheno.stats[r,grep("\\|t.stat\\|", colnames(pheno.stats))]))
							sig.edges <- rbind(sig.edges, table.row)
							}
						}
					} #end adding phenotype edges

			}
	
						
			
			
		edgelist <- matrix(c(as.vector(sig.edges[,1]), as.vector(sig.edges[,2])), ncol = 2, byrow = FALSE)
				
	
	
		net <- graph.edgelist(edgelist)
		V(net)$color <- rep("lightblue", vcount(net))
	
		#find the locations of the phenotype nodes
		pheno.locale <- match(phenotypes, V(net)$name)
		#color the phenotype nodes a different color
		if(length(pheno.locale) > 0){
			V(net)$color[pheno.locale] <- "violet"
			}


		if(class(layout.matrix) != "matrix" && !is.null(layout.matrix)){
			if(layout.matrix == "manual"){
				tkp.id <- tkplot(net)
				done <- readline(prompt = "Press return when ready:\n")
				layout.matrix <- tkplot.getcoords(tkp.id)
				tkplot.close(tkp.id)
				}
			}

		if(!is.null(layout.matrix) & length(layout.matrix[,1]) != vcount(net)){
			stop("The layout matrix does not have the same number of nodes as the network.")
			}
		
		
		if(is.null(layout.matrix)){
			
			# coord.matrix <- layout.kamada.kawai(net)
			coord.matrix <- layout.auto(net)
				
			if(separate.pheno){

				#If there are phenotype nodes, pull them out to the
				#bottom to clarify the graph a little bit
				if(length(pheno.locale) > 0){
				
					marker.locale <- c(1:vcount(net))[-pheno.locale]
		
					marker.subnet <- induced.subgraph(net, marker.locale, impl = "auto")
					# plot(marker.subnet, vertex.label = V(marker.subnet)$name)
					marker.layout <- layout.kamada.kawai(marker.subnet)

					marker.layout[,2] <- marker.layout[,2] + abs(min(marker.layout[,2]))
					coord.matrix[marker.locale,] <- marker.layout

					pheno.x.coord <- pretty(c(min(coord.matrix[,1]):max(coord.matrix[,1])), n = length(pheno.locale))[1:length(pheno.locale)]
					pheno.y.coord <- rep(min(coord.matrix[,2])*1.25, length(pheno.locale))
		
					coord.matrix[pheno.locale,1] <- pheno.x.coord
					coord.matrix[pheno.locale,2] <- pheno.y.coord
					}
				}
			}else{
			coord.matrix <- layout.matrix
			}
			
		edge.weights <- as.numeric(as.vector(sig.edges[,"Effect"]))
		#scale the edge weights and make sure they're all positive
		if(weight.edges){
			scaled.weights <- (abs(edge.weights)/max(abs(edge.weights)))*4
			}else{
			scaled.weights <- rep(1, ecount(net))
			}
		
		E(net)$weight <- scaled.weights
	
		if(!is.null(highlight.node.edges)){
			vertex.locale <- match(highlight.node.edges, V(net)$name)
			V(net)$color[vertex.locale] <- "yellow"
			if(length(which(!is.na(vertex.locale))) < length(highlight.node.edges)){
				message("I couldn't find the following nodes to highlight:")
				cat(setdiff(highlight.node.edges, V(net)$name), sep = "\n")
				}else{
				nearest.neighbors <- get.adjedgelist(net)
				edges.to.highlight <- NULL
				for(i in highlight.node.edges){
					edges.to.highlight <- unique(c(edges.to.highlight, nearest.neighbors[[i]]))
					}
					not.highlighted <- setdiff(E(net), edges.to.highlight)
					E(net)$weight[not.highlighted] <- 0.5
					E(net)$weight[edges.to.highlight] <- 8
				}

			}
	
		neg.edges <- which(edge.weights < 0)
		pos.edges <- which(edge.weights > 0)
		edge.col <- rep("gray", ecount(net))
	
		if(length(neg.edges) > 0){
			edge.col[neg.edges] <- "red"
			}
		if(length(pos.edges) > 0){
			edge.col[pos.edges] <- "green"
			}
	
		E(net)$color <- edge.col
		
		plot(net, vertex.label = V(net)$name, vertex.label.dist = 0.2, layout = coord.matrix, vertex.size = 3, edge.width = E(net)$weight)
		legend("topright", legend = c("enhancement", "repression"), col = c("green", "red"), lty = 1, lwd = 3)
		# legend("topright", legend = c("+", "-"), col = c("green", "red"), lty = 1, lwd = 3)

		results <- list(net, coord.matrix)
		names(results) <- c("net", "layout.matrix")
		invisible(results)

	
	}
