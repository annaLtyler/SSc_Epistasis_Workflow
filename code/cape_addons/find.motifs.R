#This function finds specific motifs in the gene interaction network
#if sig.motifs is set to FALSE, this finds the motifs with the smallest
#interaction coefficients. This is for comparing markers from significant
#motifs to markers from additive interactions

find.motifs <- function(data.obj, collapsed.net = TRUE, include.covar = FALSE){


	if(collapsed.net){
		total.net <- data.obj$collapsed.net
		}else{
		total.net <- data.obj$full.net
		}

		# total.net[which(total.net == 0)] <- NA
		# imageWithText(total.net, split.at.vals = TRUE, show.text = FALSE)
		# imageWithTextColorbar(total.net, split.at.vals = TRUE)

		gene.ind <- 1:min(dim(total.net))
		full.ind <- 1:max(dim(total.net))
	
	if(!include.covar){
		covar.info <- get.covar(data.obj)
		covar.names <- covar.info$covar.names
		covar.locale <- match(covar.names, colnames(total.net))
		if(length(covar.locale) > 0 && !is.na(covar.locale)){
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
	
	#find all the edges
	gene.edge.locale <- which(gene.net != 0, arr.ind = TRUE)
	
	if(length(gene.edge.locale) == 0){stop("There are no edges at this p value.")}
	
	#to find the complete three-node motifs, look at each
	#edge and verify that each gene has a connection to the
	#phenotype as well. Each phenotype will get its own set
	#of motifs
	all.triplets <- vector(mode = "list", length = num.pheno)
	names(all.triplets) <- colnames(pheno.net)
	
	all.triplet.signs <- vector(mode = "list", length = num.pheno)
	names(all.triplet.signs) <- colnames(pheno.net)
	
	for(ph in 1:num.pheno){
		triplet.list <- matrix(NA, ncol = 3, nrow = dim(gene.edge.locale)[1])
		triplet.signs <- matrix(NA, ncol = 3, nrow = dim(gene.edge.locale)[1])
		colnames(triplet.list) <- c("source", "target", "pheno")
		colnames(triplet.signs) <- c("source-target", "source-Ph", "target-Ph")
		for(g in 1:nrow(gene.edge.locale)){
			#for each edge, look for an entry in the phenotype
			#matrix corresponding to each gene in the interaction
			#to find open triplets, just change && to ||
			if(pheno.net[gene.edge.locale[g,1],ph] != 0 && pheno.net[gene.edge.locale[g,2],ph] != 0){
				triplet.list[g,] <- c(colnames(gene.net)[gene.edge.locale[g,]], colnames(pheno.net)[ph])
				triplet.signs[g,] <- c(sign(gene.net[gene.edge.locale[g,1], gene.edge.locale[g,2]]), sign(pheno.net[gene.edge.locale[g,1],ph]), sign(pheno.net[gene.edge.locale[g,2],ph]))
				}
			}
		#take out the NA rows, and the non-interactions 
		#for each phenotype
		non.na.locale <- which(!is.na(triplet.list[,1]))
		triplet.list <- triplet.list[non.na.locale,,drop=FALSE]
		triplet.signs <- triplet.signs[non.na.locale,,drop=FALSE]

		all.triplets[[ph]] <- triplet.list	
		all.triplet.signs[[ph]] <- triplet.signs
	}
	
	result <- list(all.triplets, all.triplet.signs, collapsed.net)
	return(result)

	#make more scripts for processing the motifs	
	#next step, condense A -> B and A <- B to A <-> B
	#and sort into different classes: enhancing, suppressive, 
	#main effects the same, main effects different


		
	
	
}