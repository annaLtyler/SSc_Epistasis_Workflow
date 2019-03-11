#This function counts the number of individuals
#in each pairwise genotype to help in selecting 
#pairs of markers for the pairscan


count.geno.pairs <- function(data.obj, geno.obj = NULL, genotypes = c(0, 0.5, 1)){
	
	geno <- get.geno(data.obj, geno.obj)
	num.types <- apply(geno, 2, function(x) length(unique(x)))
	too.many <- which(num.types > 3)
	too.many.mat <- geno[,too.many,drop=FALSE]
	
	binned.geno <- apply(too.many.mat, 2, function(x) bin.vector(x, genotypes))
	new.geno <- geno
	new.geno[,too.many] <- binned.geno
	
	num.pairs <- choose(ncol(new.geno), 2)
	
	pair.nums <- matrix(NA, nrow = num.pairs, ncol = 3)
	pair.counter <- 1
	for(i in 1:ncol(new.geno)){
		for(j in 2:ncol(new.geno)){
			if(i != j){
				report.progress(pair.counter, num.pairs)
				geno.counts <- table(new.geno[,i], new.geno[,j])
				pair.nums[pair.counter,] <- c(i,j, min(geno.counts))
				pair.counter <- pair.counter + 1
				}
			}
		}
		
	return(pair.nums)
}