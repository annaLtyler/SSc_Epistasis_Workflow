#This function is the same as get.network except that
#the network can be made up of the N top interaction 
#effects AND the N top main effects, or a particular 
#edge list that includes main effect. top.N.inter 
#is an integer, top.N.main is named a vector with an 
#entry for each phenotype specifying the number of 
#effects to use for each phenotype. The edge list
#should be a matrix with two columns and node names
#as the entries. The matrix must have ordered pairs
#with Source markers (not block names) in the first 
#column and Target markers in the second column. 
#One of p.or.q, top.N, or edge.list must be used.
#All parameters for selecting edges act on uncollapsed
#blocks. Collapsing into blocks may decrease the number
#of elements in the network from the top.N arguments.


get.network2 <- function(data.obj, p.or.q = NULL, top.N.inter = NULL, top.N.main = NULL, edge.list = NULL, min.std.effect = 0, collapse.linked.markers = TRUE, r.thresh = 0.5, verbose = FALSE, plot.linkage.blocks = FALSE){
	
	
	linkage.method = "genotype"

	any.net.type <- c(!is.null(p.or.q), !is.null(top.N.inter), !is.null(edge.list))
	if(sum(any.net.type) != 1){
		stop("One and only one of p.or.q, top.N, or edge.list must be specified")
		}
	
	top.N.check <- c(!is.null(top.N.inter), !is.null(top.N.main))
	if(sum(top.N.check) == 1){
		stop("Both top.N.inter and top.N.main must be specified.")
		}
	
	#get the linkage blocks based on the significant markers
	# data.obj <- linkage.blocks(data.obj, collapse.linked.markers = collapse.linked.markers, threshold.power = threshold.power, verbose = verbose, plot.results = plot.linkage.blocks)
	# method.check <- grep("eff", linkage.method)
	# if(length(method.check) > 0){linkage.method <- "effects"}
	# if(linkage.method == "effects"){
		# data.obj <- linkage.blocks.cormat(data.obj, p.or.q = p.or.q, effect.drop = effect.drop, collapse.linked.markers = collapse.linked.markers, threshold.power = threshold.power, verbose = verbose, plot.results = plot.linkage.blocks)
		# }
	# if(linkage.method == "genotype"){
		data.obj.col <- linkage.blocks.stepwise(data.obj, collapse.linked.markers = TRUE, r.thresh = r.thresh, plot.blocks = plot.linkage.blocks)
		data.obj.uncol <- linkage.blocks.stepwise(data.obj, collapse.linked.markers = FALSE, r.thresh = r.thresh, plot.blocks = plot.linkage.blocks)
		# }
	# if(linkage.method == "prominence"){
		# data.obj <- linkage.blocks.prominence(data.obj, p.or.q = p.or.q, collapse.linked.markers = collapse.linked.markers, threshold.power = threshold.power, verbose = verbose, plot.results = plot.linkage.blocks)
		# }
	
	
	if(collapse.linked.markers){
		data.obj <- data.obj.col
		blocks <- data.obj.col$linkage.blocks.collapsed
		}else{
		data.obj <- data.obj.uncol
		blocks <- data.obj.uncol$linkage.blocks.full	
		}

	full.blocks <- data.obj.uncol$linkage.blocks.full
	
	block.names <- names(blocks)	
	
	#get the data for the network
	all.net.data <- data.obj$var.to.var.p.val
	pheno.tables <- data.obj$max.var.to.pheno.influence
	phenotypes <- names(pheno.tables)	

	
	#filter the edges by whatever means is set by the user
	if(!is.null(p.or.q)){
		var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(all.net.data))))
		#filter to just the significant interactions
		net.data <- all.net.data[which(as.numeric(all.net.data[,var.sig.col]) <= p.or.q),,drop = FALSE]	
		}
	if(!is.null(top.N.inter)){
		net.data <- all.net.data[1:top.N.inter,,drop=FALSE]
		}	
	if(!is.null(edge.list)){
		#first remove the edges with phenotypes
		pheno.locale <- which(edge.list[,2] %in% phenotypes)
		not.pheno.locale <- which(!(edge.list[,2] %in% phenotypes))
		just.pheno <- edge.list[pheno.locale,]
		just.inter <- edge.list[not.pheno.locale,]
		
				
		get.edge.locale <- function(edge.pair){
			source.block.locale <- which(names(full.blocks) %in% edge.pair[1])
			target.block.locale <- which(names(full.blocks) %in% edge.pair[2])
			source.locale <- which(all.net.data[,1] %in% full.blocks[[source.block.locale]])
			target.locale <- which(all.net.data[,2] %in% full.blocks[[target.block.locale]])
			return(intersect(source.locale, target.locale))
			}
		
		if(verbose){cat("Matching edges from original network...\n")}
		edge.locale <- apply(just.inter, 1, get.edge.locale)
		net.data <- all.net.data[edge.locale,]
		}
		
	#and filter by effect size if appropriate
	net.data <- net.data[which(as.numeric(net.data[,5]) >= min.std.effect),,drop = FALSE]
	
	adj.mat <- matrix(0, ncol = length(blocks), nrow = length(blocks))
	colnames(adj.mat) <- rownames(adj.mat) <- names(blocks)
	
	block.markers <- NULL
	for(i in 1:length(blocks)){
		block.markers <- rbind(block.markers, cbind(rep(names(blocks)[i], length(blocks[[i]])), blocks[[i]]))	
		}

	
	#get the block pairs with significant effects
	get.sig.block.pair <- function(marker.pair){
		block1 <- which(block.markers[,2] == marker.pair[1])
		block2 <- which(block.markers[,2] == marker.pair[2])
		if(length(block1) > 0 && length(block2) > 0){
			return(c(block.markers[block1,1], block.markers[block2,1]))
			}else{
				return(c(0,0))
				}
		}

	sig.block.pairs <- unique(t(apply(matrix(net.data[,1:2], ncol = 2), 1, get.sig.block.pair)))

	if(length(sig.block.pairs) > 0){
	
		sig.block.pairs <- sig.block.pairs[which(sig.block.pairs[,1] != 0),,drop=FALSE]
	
		#for each pair of blocks
		get.adj.weight <- function(block.pair){
			#get all the markers in the two blocks
			all.markers1 <- blocks[[block.pair[1]]]
			all.markers2 <- blocks[[block.pair[2]]]
			
			#find the maximum weight between markers in the blocks
			#in both directions
			block1.source.locale <- which(net.data[,"Source"] %in% all.markers1)
			block2.target.locale <- which(net.data[,"Target"] %in% all.markers2)
			block1.to.block2 <- intersect(block1.source.locale, block2.target.locale)
			
			if(length(block1.to.block2) > 0){
				all.effects <- as.numeric(net.data[block1.to.block2,"Effect"])/as.numeric(net.data[block1.to.block2,"SE"])
				adj.mat[block.pair[1], block.pair[2]] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
				}
			
	
			block2.source.locale <- which(net.data[,"Source"] %in% all.markers2)
			block1.target.locale <- which(net.data[,"Target"] %in% all.markers1)
			block2.to.block1 <- intersect(block2.source.locale, block1.target.locale)
			
			if(length(block2.to.block1) > 0){
				all.effects <- as.numeric(net.data[block2.to.block1,"Effect"])/as.numeric(net.data[block2.to.block1,"SE"])
				adj.mat[block.pair[2], block.pair[1]] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
				}
				
			return(adj.mat)
			}
		
	
	
		for(i in 1:length(sig.block.pairs[,1])){
			adj.mat <- get.adj.weight(sig.block.pairs[i,])
			}
		}	
		
	#============================================================
	#internal functions for main effect placement
	#============================================================
	get.marker.block <- function(marker.name, blocks){
		return(names(unlist(lapply(blocks, function(x) which(x == marker.name)))))
		}

	place.block.inf <- function(sig.list, pheno.num, blocks){
		#find the block associated with each marker
		marker.blocks <- apply(matrix(sig.list[,1], ncol = 1), 1, get.marker.block, blocks)
		u_blocks = unique(marker.blocks)
		for(b in 1:length(u_blocks)){
			block.locale <- which(marker.blocks == u_blocks[b])
			all.effects <- as.numeric(sig.list[block.locale,"coef"])/as.numeric(sig.list[block.locale,"se"])
			pheno.mat[u_blocks[b],pheno.num] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
			}
		return(pheno.mat)
		}
	#============================================================

	#Now add the phenotypic effects continuing to use the maximum significant effect from each block
	pheno.mat <- matrix(0, nrow = length(blocks), ncol = length(phenotypes))
	colnames(pheno.mat) <- phenotypes
	rownames(pheno.mat) <- names(blocks)
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.tables[[1]]))))

	for(ph in 1:length(phenotypes)){
		if(!is.null(p.or.q)){
			sig.inf <- pheno.tables[[ph]][which(pheno.tables[[ph]][,pheno.sig.col] <= p.or.q),,drop = FALSE]
			}
		if(!is.null(top.N.main)){
			pheno.locale <- which(names(top.N.main) == phenotypes[ph])
			sig.inf <- pheno.tables[[ph]][1:top.N.main[pheno.locale],,drop = FALSE]
			}
		if(!is.null(edge.list)){
			target.locale <- which(just.pheno[,2] == phenotypes[ph])
			source.list <- just.pheno[target.locale,1]
			source.markers <- apply(matrix(source.list, ncol = 1), 1, function(x) full.blocks[[which(names(full.blocks) == x)]])
			source.locale <- which(pheno.tables[[ph]][,1] %in% source.markers)
			sig.inf <- pheno.tables[[ph]][source.locale,,drop=FALSE]
			}
		pheno.mat <- place.block.inf(sig.list = sig.inf, pheno.num = ph, blocks)
		}
	
	final.mat <- cbind(adj.mat, pheno.mat)

	if(collapse.linked.markers){
		data.obj$collapsed.net <- final.mat
		}else{
		data.obj$full.net <- final.mat	
		}
	
	return(data.obj)
	
	}






