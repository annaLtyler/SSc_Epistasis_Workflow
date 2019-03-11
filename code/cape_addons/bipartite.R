#This script generates a bipartite graph given a list of source
#nodes, a list of target nodes, and the graph
#edge.list should be an edgelist, and the nodes should be
#vectors

bipartite <- function(edge.list, source.nodes, target.nodes){

	require(igraph)	
	
	find.links <- function(edge.list, test.node, target.list){
		bi.graph <- NULL
		test.locale <- which(edge.list[,1] == test.node)
		target.locale <- which(target.list %in% edge.list[test.locale,2])
		if(length(target.locale) > 0){
			bi.graph <- rbind(bigraph, edge.list[target.locale,])
			}
		test.locale <- which(edge.list[,2] == test.node)
		target.locale <- which(target.list %in% edge.list[test.locale,1])
		if(length(target.locale) > 0){
			bi.graph <- rbind(bigraph, edge.list[target.locale,])
			}
		return(bi.graph)
		}
		
	bi.list.source <- apply(matrix(source.nodes, ncol = 1), 1, function(x) find.links(edge.list, x, target.nodes))
	bi.list.target <- apply(matrix(target.nodes, ncol = 1), 1, function(x) find.links(edge.list, x, source.nodes))
	
}