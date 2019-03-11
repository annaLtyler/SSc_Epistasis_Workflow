#This function finds motifs in the gene interaction network
#that are not epistatic. It finds the n pairs of markers 
#with significant main effects and the lowest interaction
#coefficients.
#max.effect refers to the maximum standardized effect
#size allowed for the interaction
#min.pval refers to the minimum p value allowed for the
#interaction.
#one or both may be specified
#The p value for the main effects is determined by the 
#input network.
#right now this only works for the full network
#all names should be in terms of chromosomal regions
#as on the network objects

find.add.motifs <- function(data.obj, max.effect = 0.5, min.pval = 0.5, include.covar = FALSE){


	if(collapsed.net){
		total.net <- data.obj$collapsed.net
		name.key <- data.obj$linkage.blocks.collapsed
		}else{
		total.net <- data.obj$full.net
		name.key <- data.obj$linkage.blocks.full
		}

		gene.ind <- 1:min(dim(total.net))
		full.ind <- 1:max(dim(total.net))
	
	if(!include.covar){
		covar.ind <- unique(as.vector(apply(data.obj$covar.for.pairscan, 2, function(x) which(x == 1))))
		covar.names <- data.obj$marker.names[match(covar.ind, colnames(data.obj$geno))]
		covar.locale <- match(covar.names, colnames(total.net))
		if(length(covar.locale) > 0){
			gene.ind <- gene.ind[-covar.locale]
			full.ind <- full.ind[-covar.locale]
			}
		}
		
	pheno.names <- colnames(data.obj$pheno)
	pheno.ind <- match(pheno.names, colnames(total.net)) 
	
	#pull out the gene network
	gene.net <- total.net[gene.ind, gene.ind]
	#and the edges to the phenotypes
	pheno.net <- total.net[gene.ind, pheno.ind]
	num.pheno <- dim(pheno.net)[2]
	
	#find all pairs of markers with significant main effects
	#in any of the phenotypes
	
	non.sig.interactions <- data.obj$var.to.var.p.val
	#Find all the interaction terms with an effect size less than or equal to the user input
	if(!is.null(max.effect)){
		non.sig.interactions <- non.sig.interactions[which(non.sig.interactions[,5] <= max.effect),]
		}
	if(!is.null(min.pval)){
		non.sig.interactions <- non.sig.interactions[which(non.sig.interactions[,7] >= min.pval),]		
		}
	

	both.ways <- function(marker.pair){
		way.one <- intersect(which(non.sig.interactions[,1] == marker.pair[1]), which(non.sig.interactions[,2] == marker.pair[2]))
		way.two <- intersect(which(non.sig.interactions[,1] == marker.pair[2]), which(non.sig.interactions[,2] == marker.pair[1]))
		if(length(way.one) > 0 & length(way.two) > 0){
			return(TRUE)
			}else{
			return(FALSE)
			}
		}
	
	both.non.sig <- apply(non.sig.interactions[,1:2], 1, both.ways)
	
	if(length(which(both.non.sig)) == 0){
		stop("There are no motifs matching these criteria. Please loosen the restrictions.\n")
		}
	both.non.sig.with.main <- non.sig.interactions[which(both.non.sig),]
	
	#filter out covariates from this table
	marker.names <- apply(both.non.sig.with.main[,1:2], 2, function(x) data.obj$marker.names[match(x, colnames(data.obj$geno))])
	if(!include.covar){
		covar.locale <- unique(as.vector(apply(marker.names, 2, function(x) which(x %in% covar.names))))
		if(length(covar.locale) > 0){
			both.non.sig.with.main <- both.non.sig.with.main[-covar.locale,]
			}
		}
	
	#to find the complete three-node motifs, look at each
	#edge and verify that each gene has a connection to the
	#phenotype as well. Each phenotype will get its own set
	#of motifs
	all.triplets <- vector(mode = "list", length = num.pheno)
	names(all.triplets) <- colnames(pheno.net)
	
	all.triplet.signs <- vector(mode = "list", length = num.pheno)
	names(all.triplet.signs) <- colnames(pheno.net)
	
	for(ph in 1:num.pheno){
		triplet.list <- matrix(NA, ncol = 3, nrow = dim(both.non.sig.with.main)[1])
		triplet.signs <- matrix(NA, ncol = 3, nrow = dim(both.non.sig.with.main)[1])
		colnames(triplet.list) <- c("source", "target", "pheno")
		colnames(triplet.signs) <- c("source-target", "source-Ph", "target-Ph")
		for(g in 1:dim(both.non.sig.with.main)[1]){
			#find the main effects for each pair using the input network
			marker1 <- both.non.sig.with.main[g,1]
			marker1.block <- names(unlist(lapply(name.key, function(x) which(x == as.vector(marker1)))))
			
			marker2 <- both.non.sig.with.main[g,2]
			marker2.block <- names(unlist(lapply(name.key, function(x) which(x == as.vector(marker2)))))
			
			marker1.locale <- which(rownames(pheno.net) == marker1.block)
			marker2.locale <- which(rownames(pheno.net) == marker2.block)
			#if both markers have significant main effects
			if(pheno.net[marker1.locale,ph] != 0 && pheno.net[marker2.locale,ph] != 0){
				#record the marker and phenotype names
				triplet.list[g,] <- c(marker1.block, marker2.block, colnames(data.obj$pheno)[ph])
				#and the main effect signs
				triplet.signs[g,] <- c(0, sign(pheno.net[marker1.locale,ph]), sign(pheno.net[marker2.locale,ph]))
				}
			}
		#take out the NA rows, and the non-interactions 
		#for each phenotype
		non.na.locale <- which(!is.na(triplet.list[,1]))
		triplet.list <- triplet.list[non.na.locale,]
		triplet.signs <- triplet.signs[non.na.locale,]

		all.triplets[[ph]] <- triplet.list	
		all.triplet.signs[[ph]] <- triplet.signs
	}
	
	result <- list(all.triplets, all.triplet.signs, collapsed.net)
	return(result)
	
}