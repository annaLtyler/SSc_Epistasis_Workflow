#This function takes in two data.objects: data.obj1 has 
#a desired network gotten through get.network with a p.or.q 
#value. data.obj2 is another data object with its own main
#effects and interaction values but whose network you want
#to match to that in data.obj1. 
#The function generates a network for data.obj2 with 
#parameters gathered from the data.obj1 network.This 
#function can match the network either through getting 
#the same number of edges as in the original, or by using 
#the exact edges as in the original. This would be used if 
#you are trying to compare a cape network with an epistasis 
#network and you want as similar networks as possible

match.net <- function(data.obj1, data.obj2, orig.pval = 0.01, match.type = c("N", "edges"), verbose = FALSE){
	
	data.obj1 <- get.network(data.obj1, p.or.q = orig.pval, collapse.linked.markers = TRUE)
	data.obj1 <- get.network(data.obj1, p.or.q = orig.pval, collapse.linked.markers = FALSE)
	
	full.net <- data.obj1$full.net
	num.markers <- min(dim(full.net))
	
	marker.locale <- 1:num.markers
	pheno.locale <- (num.markers+1):max(dim(full.net))
	
	just.inter <- full.net[,marker.locale]
	just.main <- full.net[,pheno.locale]
	
	if(match.type == "N"){
		n.inter <- length(which(just.inter != 0))
		n.main <- apply(just.main, 2, function(x) length(which(x != 0)))
		
		data.obj2 <- get.network2(data.obj2, top.N.inter = n.inter, top.N.main = n.main, collapse.linked.markers = TRUE, verbose = verbose)
		data.obj2 <- get.network2(data.obj2, top.N.inter = n.inter, top.N.main = n.main, collapse.linked.markers = FALSE, verbose = verbose)
		}
	
	if(match.type == "edges"){
		all.edges <- which(full.net != 0, arr.ind = TRUE)
		edge.mat <- matrix(NA, nrow = nrow(all.edges), ncol = ncol(all.edges))
		edge.mat[,1] <- rownames(full.net)[all.edges[,1]]
		edge.mat[,2] <- colnames(full.net)[all.edges[,2]]
		
		if(verbose){cat("Matching Full Network...\n")}
		data.obj2 <- get.network2(data.obj = data.obj2, edge.list = edge.mat, collapse.linked.markers = TRUE, verbose = verbose)
		if(verbose){cat("Matching Collapsed Network...\n")}
		data.obj2 <- get.network2(data.obj2, edge.list = edge.mat, collapse.linked.markers = FALSE, verbose = verbose)
		}
	
	return(data.obj2)
	
}