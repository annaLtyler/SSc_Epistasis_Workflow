#this function determines whether A given marker pair
#had a standard kinship correction or a kinship correction
#with forced positive covariance

get.correction.status <- function(data.obj, kinship.obj, pheno.num, pair){
	
	kin <- kinship.obj[[pheno.num]]

	chr.pairs <- strsplit(names(kinship.obj[[1]]), ",")
	chr.pair.mat <- cbind(unlist(lapply(chr.pairs, function(x) x[1])), unlist(lapply(chr.pairs, function(x) x[2])))

	rejected.markers <- data.obj$rejected.markers

	marker.names <- colnames(data.obj$geno)
	if(length(rejected.markers) > 0){
		rej.marker.locale <- which(marker.names %in% rejected.markers)
		marker.names <- marker.names[-rejected.markers]
		}

	marker.locale <- match(pair, marker.names)
	marker.chr <- matrix(data.obj$chromosome[marker.locale], nrow = 1)
	chr.locale <- intersect(which(chr.pair.mat[,1] == marker.chr[1]), which(chr.pair.mat[,2] == marker.chr[2]))
				
	kindat <- kin[[chr.locale]]
	return(kindat$correction.flag)
	
}