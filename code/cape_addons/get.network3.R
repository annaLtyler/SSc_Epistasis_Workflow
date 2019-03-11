#This function collapses a network based on
#the full network. Nodes with the same connections
#and same main effects, they are combined into
#a single node


get.network3 <- function(data.obj){
	
	
	full.net <- data.obj$full.net
	num.pheno <- dim(data.obj$pheno)[2]
	
	#go through each row of the network
	#and find nodes with the same sign 
	#main effects and interactions
	#whittle away the full network matrix
	#as we go
	
	signed.net <- matrix(0, ncol = dim(full.net)[2], nrow = dim(full.net)[1])
	signed.net[which(full.net > 0)] <- 1
	signed.net[which(full.net < 0)] <- -1
	rownames(signed.net) <- rownames(full.net)
	colnames(signed.net) <- colnames(full.net)
	
	marker.names <- data.obj$marker.names[which(colnames(cross$geno) %in% colnames(full.net))]
	
	#for each node, find its ingoing and outgoing connections
	get.in.out <- function(node.num){
		incoming <- signed.net[,node.num]; names(incoming) <- rep("in", length(incoming))
		outgoing <- signed.net[node.num,]; names(outgoing) <- rep("out", length(outgoing))
		in.out <- c(incoming, outgoing)
		return(in.out)
		}	
	
	inout.mat <- apply(matrix(1:dim(signed.net)[1], ncol = 1), 1, get.in.out)
	colnames(inout.mat) <- marker.names
	
	#find the unique columns of the inout matrix
	unique.nodes <- unique(inout.mat, MARGIN = 2)

	#rebuild a new compressed matrix
	unique.node.locale <- which(marker.names %in% colnames(unique.nodes))
	compressed.mat <- full.net[unique.node.locale,unique.node.locale]
	
	pheno.locale <- which(colnames(full.net) %in% colnames(data.obj$pheno))
	compressed.mat <- cbind(compressed.mat, full.net[unique.node.locale,pheno.locale])
		
	#go through the unique nodes, find the nodes they 
	#share effects with so we can name each unique node
	
	
	for(i in 1:dim(unique.nodes)[2]){
		same.nodes <- names(which(apply(inout.mat, 2, function(x) identical(unique.nodes[,i], x))))
		split.name <- strsplit(same.nodes, "_")
		new.name.chr <- paste(unique(sapply(split.name, function(x) x[1])), collapse = ",")
		new.name.id <- paste(sapply(split.name, function(x) x[2]), collapse = ",")
		new.name <- paste(c(new.name.chr, new.name.id), collapse = "_")
		colnames(compressed.mat)[i] <- new.name
		rownames(compressed.mat)[i] <- new.name
		}
	
	data.obj$collapsed.net <- compressed.mat
	return(data.obj)	
	}

