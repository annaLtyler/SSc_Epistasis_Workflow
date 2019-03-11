#This function condenses the genome into the first PC of each 
#linkage block. I am currently testing whether this is a 
#reasonable way to reduce the amount of testing, although it
#might reduce the signal too much

compress.genome <- function(data.obj){
	
	not.covar <- which(data.obj$chromosome != 0)
	just.geno <- data.obj$geno[,not.covar]
	is.na <- which(is.na(just.geno))
	if(length(is.na) > 0){stop("There cannot be missing data in the genotypes. Use impute.geno() to impute genotypes.")}

	geno <- apply(data.obj$geno, 2, function(x) x - mean(x, na.rm = TRUE))

	linkage.blocks <- data.obj$linkage.blocks.collapsed
	
	#get the covariate names. If we are looking at a covariate,
	#just put it straight in.
	covar.names <- data.obj$marker.names[which(data.obj$chromosome == 0)]
	new.geno <- matrix(NA, nrow = nrow(geno), ncol = length(linkage.blocks))


	marker.location <- rep(NA, length(linkage.blocks))
	marker.chromosome <- rep(NA, length(linkage.blocks))
	for(i in 1:length(linkage.blocks)){
		is.covar <- which(covar.names == names(linkage.blocks)[i])
		if(length(is.covar) > 0){
			covar.locale <- which(data.obj$marker.names == names(linkage.blocks)[i])
			new.geno[,i] <- geno[,covar.locale]
			marker.location[i] <- data.obj$marker.location[covar.locale]
			marker.chromosome[i] <- data.obj$chromosome[covar.locale]
			}else{
			geno.locale <- match(linkage.blocks[[i]], colnames(geno))
			marker.location[i] <- mean(data.obj$marker.location[geno.locale])
			marker.chromosome[i] <- unique(data.obj$chromosome[geno.locale])
			geno.block <- geno[,geno.locale,drop=FALSE]
			new.geno[,i] <- svd(geno.block)$u[,1]
			}
		}
	
	colnames(new.geno) <- 1:dim(new.geno)[2]
	data.obj$geno <- new.geno
	data.obj$marker.names <- names(linkage.blocks)
	data.obj$chromosome <- marker.chromosome
	data.obj$marker.location <- marker.location
	
	return(data.obj)
	
}