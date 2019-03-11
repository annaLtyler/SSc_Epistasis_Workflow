#This script converts R/qtl objects to the cape data object

rqtl2cape <- function(rqtl.object, imputed.markers = TRUE){
	
	pheno <- as.matrix(rqtl.object$pheno)
	rownames(pheno) <- 1:dim(pheno)[1]


	geno.obj <- rqtl.object$geno
	
	geno.obj$X <- NULL
	
	num.ind <- dim(geno.obj[[1]][[1]])[1]
	
	#check to see if there are imputed psuedomarkers
	is.imputed <- !is.null(rqtl.object$geno[[1]]$prob)
	
	#If not, just use the genotype
	if(!imputed.markers){
		num.markers <- sum(sapply(geno.obj, function(x) dim(x[[1]])[2]))

		geno.mat <- matrix(NA, nrow = num.ind, ncol = num.markers)	
		marker.names <- rep(NA, num.markers)
		chromosome <- rep(NA, num.markers)
		marker.pos <- rep(NA, num.markers)
	
		for(i in 1:length(geno.obj)){
			first.na <- min(which(apply(geno.mat, 2, function(x) length(which(is.na(x)))) == dim(geno.mat)[1]))
			chr.geno <- as.matrix(geno.obj[[i]][[1]])
			num.chr.markers <- dim(chr.geno)[2]
			geno.mat[,first.na:(first.na+num.chr.markers-1)] <- chr.geno
			
			marker.names[first.na:(first.na+num.chr.markers-1)] <- colnames(chr.geno)
			marker.pos[first.na:(first.na+num.chr.markers-1)] <- as.vector(geno.obj[[i]][[2]])
			chromosome[first.na:(first.na+num.chr.markers-1)] <- rep(names(rqtl.object$geno)[i], num.chr.markers)
			}
		colnames(geno.mat) <- 1:dim(geno.mat)[2]
		
		
		#make sure the genotypes are between 0 and 1
		u_geno <- sort(unique(geno.mat[which(!is.na(geno.mat))]))
		new.geno <- seq(0,1,(1/(length(u_geno)-1)))
		
		for(i in 1:length(u_geno)){
			geno.mat[which(geno.mat == u_geno[i])] <- new.geno[i]
			}
			
		}else{
			if(!is.imputed){
				stop("I can't find any imputed markers from calc.genoprob in R/qtl.")
				}
				
			#otherwise we are dealing with pseudomarkers, which is a bit more complicated
			
			#get the number of markers and pseudomarkers
			num.markers <- sum(sapply(geno.obj, function(x) dim(x$prob)[2]))
			
			geno.mat <- matrix(NA, nrow = num.ind, ncol = num.markers)	
			marker.names <- rep(NA, num.markers)
			chromosome <- rep(NA, num.markers)
			marker.pos <- rep(NA, num.markers)

			#get the possible genotypes for each marker
			genotypes <- unique(unlist(lapply(geno.obj, function(x) dimnames(x$prob)[[3]])))
			geno.probs <- seq(0,1,(1/(length(genotypes)-1)))

			for(i in 1:length(geno.obj)){
				
				first.na <- min(which(apply(geno.mat, 2, function(x) length(which(is.na(x)))) == dim(geno.mat)[1]))
				
				pseudomarkers <- geno.obj[[i]]$prob
				num.chr.markers <- dim(pseudomarkers)[2]
				marker.names[first.na:(first.na+num.chr.markers-1)] <- paste("Chr", names(rqtl.object$geno)[i], "_", names(attr(pseudomarkers, which = "map")), sep = "")
				marker.pos[first.na:(first.na+num.chr.markers-1)] <- as.vector(attr(pseudomarkers, which = "map"))	
				chromosome[first.na:(first.na+num.chr.markers-1)] <- rep(names(rqtl.object$geno)[i], num.chr.markers)
				
				
				#for each position, find the most likely
				#genotype and use the associated probability
				
				for(j in 1:num.chr.markers){
					ind.geno <- pseudomarkers[,j,]
					called.geno <- apply(ind.geno, 1, function(x) sum(x*geno.probs))
					geno.mat[,(first.na+j-1)] <- called.geno
					}
				}
			
			
			}
	
		marker.num <- 1:dim(geno.mat)[2]
		colnames(geno.mat) <- marker.num
		rownames(geno.mat) <- rownames(pheno)
		final.object <- list(pheno, geno.mat, chromosome, marker.names, marker.num, marker.pos)
		names(final.object) <- c("pheno", "geno", "chromosome", "marker.names", "marker.num", "marker.location")
	
		return(final.object)

	
}