#This function builds a gene-gene network out of top-ranked genes 
#in the object returned from get.candidate.genes


candidate.gene.network <- function(int.loci, candidate.gene.obj, vote.min = 100){
	
	require(igraph)
	require(RColorBrewer)	
	
	cols <- brewer.pal(12, "Set3")
	
	get.top.genes <- function(snp.name){
		top.gene <- NA
		snp.locale <- which(names(candidate.gene.obj) == snp.name)
		if(length(snp.locale) > 0){
			top.gene.locale <- which(candidate.gene.obj[[snp.locale]][,3] >= 100)
			if(length(top.gene.locale) == 0){
				top.gene.locale <- 1
				}
			top.gene <- paste0(candidate.gene.obj[[snp.locale]][top.gene.locale,1], collapse = ",")
			}
		return(top.gene)
		}	
		
	get.unique.vertices <- function(vertex.names){
		#find all vertices with multiple genes
		multi.gene.locale <- grep(pattern = ",", vertex.names)
		split.genes <- strsplit(vertex.names[multi.gene.locale1], ",")
		u_genes <- unique(unlist(split.genes))
		u_gene.locale <- lapply(u_genes, function(x) grep(x, vertex.names))
		names(u_gene.locale) <- u_genes
		
		#Vertices that share genes get glommed together		
		#if any entries in u_gene.locale are shared, between
		#genes, the genes get grouped.
		gene1 <- 1
		new.vertices <- list(u_genes[gene1])
		for(i in 2:length(u_gene.locale)){
			gene2 <- gene1 + 1
			shared.vertices <- length(intersect(u_gene.locale[[gene1]], u_gene.locale[[gene2]]))
			if(length(shared.vertices) > 0){
				new.vertices[[length(new.vertices)]] <- paste(new.vertices[[length(new.vertices)]], u_genes[gene2], sep = ",")
				}else{
				new.vertices[[(length(new.vertices)+1)]] <- u_genes[gene2]
				}
			gene1 <- gene1 + 1
			}
	
		return(new.vertices)
		}
		
	replace.vertex.names <- function(vertex.names, grouped.vertex){
		split.group <- strsplit(grouped.vertex, ",")[[1]]
		name.locale <- unique(unlist(lapply(split.group, function(x) grep(x, vertex.names))))
		vertex.names[name.locale] <- grouped.vertex
		return(vertex.names)
		}
		
		
	source.genes <- unlist(lapply(int.loci[,1], get.top.genes))
	target.genes <- unlist(lapply(int.loci[,2], get.top.genes))
	
	edge.list <- cbind(source.genes, target.genes)
	na.locale <- which(is.na(edge.list), arr.ind = TRUE)
	if(nrow(na.locale) > 0){
		edge.list <- edge.list[-na.locale,]
		}

	u_edge.list <- unique(edge.list)
	
	grouped.vertices <- get.unique.vertices(vertex.names = c(u_edge.list[,1], u_edge.list[,2]))
	
	for(i in 1:length(grouped.vertices)){
		u_edge.list[,1] <- replace.vertex.names(u_edge.list[,1], grouped.vertices[[i]])
		u_edge.list[,2] <- replace.vertex.names(u_edge.list[,2], grouped.vertices[[i]])
		}
	
	#reduce to unique edges again
	u_edge.list <- unique(u_edge.list)
		
	net <- graph_from_edgelist(as.matrix(u_edge.list), directed = TRUE)

	net.entrez <- getBM(c("external_gene_name", "entrezgene"), "external_gene_name", V(net)$name, mart)
	node.expr <- lapply(net.entrez[,2], gene.expr)
	names(node.expr) <- net.entrez[,1]
	node.diff <- unlist(lapply(node.expr, function(x) x[1] - x[2]))

	diff.order <- match(V(net)$name, net.entrez[,1])
	cbind(net.entrez[diff.order,], V(net)$name)

	expr.cols <- node.diff[diff.order]
	low.vals <- which(as.numeric(expr.cols) < 0)
	high.vals <-which(as.numeric(expr.cols) > 0) 
	expr.cols[low.vals] <- "blue"
	expr.cols[high.vals] <- "yellow"
	expr.cols[which(is.na(expr.cols))] <- "gray"
	
	deg <- degree(net)

	num.com <- rep(NA, 100)
	for(i in 1:100){
	com <- cluster_walktrap(net, steps = i)$membership
	num.com[i] <- length(unique(com))
	}
	
	com <- cluster_walktrap(net, steps = 8)$membership


	pdf("SSc.Gene.Network.by.Nearest.Gene.pdf", width = 10, height = 10)
	plot(net, vertex.size = sqrt(deg)*2, vertex.color = cols[com], vertex.label.dist = 0, vertex.label.cex = 0.5, edge.arrow.size = 0.3)
	legend("topleft", fill = cols[1:max(com)], legend = paste("community", 1:max(com)), cex = 1)

	plot(net, vertex.size = sqrt(deg)*2, vertex.color = expr.cols, vertex.label.dist = 0, vertex.label.cex = 0.5, edge.arrow.size = 0.3)
	legend("topleft", fill = c("blue", "yellow", "gray"), legend = c("low in SSc", "high in SSc", "no data"), cex = 1)


	
	gene.groups <- lapply(1:max(com), function(x) V(net)$name[which(com == x)])
	for(i in 1:length(gene.groups)){
		print(i)
		group.genes <- unlist(strsplit(gene.groups[[i]], ","))
		enrichment <- gprofiler(group.genes, "hsapiens")
		plot.enrichment(enrichment, 20, plot.label = paste("Community", i, "\n", paste(group.genes, collapse = ", ")), text.size = 1)
		}

	dev.off()
}