#This function builds a network out of just the motifs
#This network can then be analyzed in igraph for 
#various statistics

build.motif.net <- function(motif.obj){

	num.pheno <- length(motif.obj[[2]])
	nets <- vector(mode = "list", length = num.pheno)
	names(nets) <- names(motif.obj[[2]])
	
	for(i in 1:num.pheno){
		nodes <- unique(c(motif.obj[[1]][[i]][,1], motif.obj[[1]][[i]][,2], motif.obj[[1]][[i]][,3]))
		net <- matrix(0, length(nodes), length(nodes))
		colnames(net) <- nodes
		rownames(net) <- nodes
		#go through each edge and add it to the network	
		for(j in 1:dim(motif.obj[[1]][[i]])[1]){
			#add interaction
			net[motif.obj[[1]][[i]][j,1], motif.obj[[1]][[i]][j,2]] <- motif.obj[[2]][[i]][j,1]
			#add source pheno
			net[motif.obj[[1]][[i]][j,1], motif.obj[[1]][[i]][j,3]] <- motif.obj[[2]][[i]][j,2]
			net[motif.obj[[1]][[i]][j,2], motif.obj[[1]][[i]][j,3]] <- motif.obj[[2]][[i]][j,3]
			}
			
		nets[[i]] <- net
		}
	
	return(nets)
	
	
}