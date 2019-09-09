#This function groups ontology terms based on the ontology structure
#obo <- readRDS("~/Documents/Data/Ontologies/GO.RData")
#terms <- c("collagen metabolic process", "extracellular matrix structural constituent", "vasoconstriction",
#"collagen trimer", "wound healing")
#terms <- c("GO:0005201", "GO:0007275", "GO:0007155", "GO:0001228", "GO:0031012", "GO:0006811", "GO:0010469", "GO:0005509", "GO:0043565", "GO:0045944", "GO:0005829", "GO:0005515", "GO:0007399", "GO:0007507", "GO:0000122", "GO:0001085", "GO:0005856", "GO:0030198", "GO:0007156", "GO:0007010", "GO:0005737", "GO:0043065", "GO:0035556", "GO:0019825", "GO:0005886", "GO:0009986", "GO:0032963", "GO:0015629", "GO:0007009", "GO:0008083", "GO:0000978", "GO:0008361", "GO:0002020", "GO:0046872", "GO:0042310", "GO:0005788", "GO:0048407", "GO:0009611", "GO:0005178", "GO:0008201", "GO:0014068", "GO:0007266", "GO:0005581", "GO:0030175", "GO:0010934", "GO:0043393", "GO:0001501", "GO:0007267", "GO:0004888", "GO:0005576", "GO:0007186")

group.onto <- function(terms, obo, cut.level = NULL, plot.results = TRUE, 
id.type = c("name", "ID"), max.cluster.size = 5, min.cluster.size = 2){

	#if we are given ID's translate to names first
	if(id.type == "ID"){
		all.id <- sapply(obo, function(x) x[1,2])	
		all.names <- sapply(obo, function(x) x[which(x[,1] == "name:")[1],2])
		term.names <- sapply(terms, function(x) all.names[which(all.id == x)[1]])
		terms <- term.names
		}

	parent.list <- lapply(terms, function(x) get.parents(x, obo, out.type = "name"))
	
	#create a table of all terms and their parents. This
	#allows us to cut at a particular level
	chain.length <- max(sapply(parent.list, length))
	parent.table <- matrix(NA, nrow = chain.length, ncol = length(parent.list))
	colnames(parent.table) <- terms
	for(i in 1:length(parent.list)){
	    	n.terms <- length(parent.list[[i]])
    		rev.terms <- rev(parent.list[[i]])
	    parent.table[1:length(rev.terms),i] <- rev.terms
	}
	
	if(!is.null(cut.level)){
		cut.table <- parent.table[cut.level:nrow(parent.table),]
	}else{
		cut.table <- parent.table
	}
	
	remove.words <- c("the", "to", "of", "in", "by", "via")
	#calculate jaccard indices between the bags of words from pairs of list elements
	split.words <- function(go.term.list){
		all.words <- unique(unlist(strsplit(go.term.list, " ")))
		all.words <- unique(unlist(strsplit(all.words, "-")))
		all.words <- setdiff(all.words, remove.words)
		all.words <- all.words[which(!is.na(all.words))]
		return(all.words)
	}
	
	split.terms <- lapply(1:ncol(cut.table), function(x) split.words(cut.table[,x]))
	non.zero <- which(sapply(split.terms, length) > 0)
	split.terms <- split.terms[non.zero]
		
	term.pairs <- pair.matrix(1:length(split.terms))
	pair.jaccard <- matrix(0, nrow = length(split.terms), length(split.terms))
	colnames(pair.jaccard) <- rownames(pair.jaccard) <- terms[non.zero]

	for(i in 1:nrow(term.pairs)){
		term1 <- term.pairs[i,1]
		term2 <- term.pairs[i,2]
		pair.jaccard[term1,term2] <- jaccard.ind(split.terms[[term1]], split.terms[[term2]])
	}
	
	if(plot.results){
		pheatmap(pair.jaccard, show_rownames = FALSE, show_colnames = FALSE)
		}
		
	net <- graph_from_adjacency_matrix(pair.jaccard, mode = "undirected", weighted = TRUE)
	clust <- iter.cluster2(net, max.mod.size = max.cluster.size, min.mod.size = min.cluster.size, 
	sep = ":")

	ordered.rows <- unlist(lapply(rownames(clust), function(x) unlist(strsplit(x, ":"))))
	return(ordered.rows)
	
}