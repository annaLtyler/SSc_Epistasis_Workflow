#This function takes in an IMP network
#and pulls out the genes in the top N
#clusters by size.
#The clusters can then be analyzed for
#GO term enrichment

get.genes.in.cluster <- function(imp.net, num.clusters = 1, plot.clusters = FALSE, v.labels = NULL, net.layout = NULL){
	
	require(igraph)
	clust <- clusters(imp.net)
	min.weight <- min(E(imp.net)$weight)
	all.cluster.sizes <- clust$csize
	sorted.cluster.sizes <- sort.int(all.cluster.sizes, decreasing = TRUE, index.return = TRUE)
	
	orig.cols <- V(imp.net)$color
	all.clust.cols <- V(imp.net)$color
	
	if(is.null(net.layout)){
		net.layout <- layout.auto(imp.net)
		}

	select.num <- sorted.cluster.sizes$ix[1:num.clusters]
	
		
	if(plot.clusters){
		pdf(paste("Clusters.thresh.", min.weight, ".pdf", sep = ""), width = 20, height = 20)
		# pdf(paste("Clusters.thresh.", min.weight, ".pdf", sep = ""), width = 8, height = 8)
		# pdf(paste("Clusters.thresh.", min.weight, ".pdf", sep = ""))
		plot(imp.net, layout = net.layout, vertex.size = 2, vertex.label = NA, main = paste("The Network"))
		plot(imp.net, layout = net.layout, vertex.size = 2, vertex.label = v.labels, main = paste("The Network"), vertex.label.dist = 0.15, vertex.label.cex = 0.8)
		}
	for(i in 1:length(select.num)){
		memb.locale <- which(clust$membership == select.num[i])
		clust.memb <- V(imp.net)$name[memb.locale]
		write.table(clust.memb, paste("Cluster.", i, ".Membership.Thresh.", min.weight, ".txt", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
		
		if(plot.clusters){
			V(imp.net)$color[memb.locale] <- "green"
			all.clust.cols[memb.locale] <- "green"
			plot(imp.net, layout = net.layout, vertex.size = 2, vertex.label = NA, main = paste("Cluster", i, sep = " "))
			V(imp.net)$color <- orig.cols
			}
		
		}
		if(plot.clusters){
			plot(imp.net, layout = net.layout, vertex.size = 2, vertex.color = all.clust.cols, vertex.label = NA, main = paste("All Clusters"))
			dev.off()
			}
	
	invisible(net.layout)
}