#This function splits the data object network
#into adjacency matrices to be used in motif.enrichment()


build.adj.mats <- function(data.obj, collapsed.net = TRUE, include.covar = FALSE, use.eigentraits = FALSE){

	if(collapsed.net){
		total.net <- data.obj$collapsed.net
		}else{
		total.net <- data.obj$full.net
		}

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
		
	if(use.eigentraits){
		pheno.names <- names(data.obj$singlescan.results)
		}else{
		pheno.names <- colnames(data.obj$pheno)	
		}
	pheno.ind <- match(pheno.names, colnames(total.net)) 
		

	#pull out the gene network
	gene.net <- total.net[gene.ind, gene.ind]
	#and the edges to the phenotypes
	pheno.net <- total.net[gene.ind, pheno.ind]
	num.pheno <- dim(pheno.net)[2]
	
	#build the edge network, this is duplicated
	#for each phenotype
	nodes <- colnames(gene.net)
	net <- matrix(0, length(nodes), length(nodes))
	colnames(net) <- nodes
	rownames(net) <- nodes
	#put in signed edges for each significant edge
	gene.net[which(gene.net < 0)] <- -1
	gene.net[which(gene.net > 0)] <- 1

	nets <- vector(mode = "list", length = num.pheno)
	names(nets) <- pheno.names
	
	for(i in 1:num.pheno){
		#add the phenotype column to the gene network
		pheno.col <- pheno.net[,i]
		pheno.col[which(pheno.col > 0)] <- 1
		pheno.col[which(pheno.col < 0)] <- -1
		net <- cbind(gene.net, pheno.col)
		net <- rbind(net, rep(0, dim(net)[2]))
		colnames(net)[dim(net)[2]] <- pheno.names[i]
		rownames(net)[dim(net)[1]] <- pheno.names[i]
		nets[[i]] <- net
		}
	
	return(nets)
	
	
}