#this function removes markers with low representation
#and individuals with to few snps
#the marker.min and ind.min are both the minimum percent
#data an individual or marker should have

remove.bad.ind.markers <- function(data.obj, marker.min = 90, ind.min = 90){
	
	geno <- data.obj$geno
	
	#count the number of markers each individual has
	perc.markers <- apply(geno, 2, function(x) length(which(!is.na(x))))/dim(geno)[1]

	
	low.markers <- which(perc.markers < (marker.min/100))

	if(length(low.markers) > 0){
		cat("These markers were removed:\n")
		cat(names(low.markers), sep = "\n")
		}
		
	markers.to.keep <- which(perc.markers >= (marker.min/100))
	new.geno <- geno[,markers.to.keep]
	new.chr <- data.obj$chromosome[markers.to.keep]
	new.names <- data.obj$marker.names[markers.to.keep]
	new.loc <- data.obj$marker.location[markers.to.keep]
	
	
	#count the number of individuals genotyped at each marker
	perc.ind <- apply(new.geno, 1, function(x) length(which(!is.na(x))))/dim(new.geno)[2]

	low.ind <- which(perc.ind < (ind.min/100))


	if(length(low.ind) > 0){
		cat("These individuals were removed:\n")
		cat(low.ind, sep = "\n")		
		}

	new.geno <- new.geno[which(perc.ind >= (ind.min/100)), ]

	
	data.obj$geno <- new.geno
	data.obj$chromosome <- new.chr
	data.obj$marker.names <- new.names
	data.obj$marker.locl <- new.loc
	
	return(data.obj)
}