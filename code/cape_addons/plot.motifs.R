#This function takes in a network motif object and
#a data object and plots just the network motifs
#for each phenotype

plot.motifs <- function(data.obj, motif.obj, filename = "Motifs.By.Phenotype.pdf"){
	
	library(igraph)
	
	phenos <- names(motif.obj[[1]])
	
	pdf(filename)
	for(ph in 1:length(phenos)){
		net.elements <- unique(c(motif.obj[[1]][[ph]][,1], motif.obj[[1]][[ph]][,2]))
	
		#change this to be flexible later
		if(motif.obj[[3]]){
			super.net <- data.obj$collapsed.net
			}else{
			super.net <- data.obj$full.net	
			}
		pheno.locale <- which(colnames(super.net) == phenos[ph])
		element.order <- match(net.elements, colnames(super.net))
		adj.mat <- super.net[element.order, c(element.order, pheno.locale)]
		#add a row for the phenotype
		adj.mat <- rbind(adj.mat, rep(0, dim(adj.mat)[2]))
		net <- graph.adjacency(adj.mat, weighted = TRUE)
		E(net)$color <- rep(NA, ecount(net))
		E(net)$color[which(E(net)$weight < 0)] <- "red"
		E(net)$color[which(E(net)$weight > 0)] <- "green"
		E(net)$width = 4
		V(net)$size <- 3
		
		# quartz();plot(net, layout = layout.circle(net))
		plot(net, layout = layout.graphopt(net))
	}
	dev.off()
	
	
}