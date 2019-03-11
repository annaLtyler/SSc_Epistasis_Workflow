#This function identifies combinations of phenotypes
#with the most individuals

max.pheno <- function(data.obj, num.traits = 3, num.samples = 100){
	
	pheno <- data.obj$pheno

	# num.missing <- apply(pheno, 2, function(x) length(which(is.na(x))))
	# pheno.order <- order(num.missing)
	# num.missing[pheno.order]
	
	total.ind <- rep(NA, num.samples)
	sampled.traits <- matrix(NA, nrow = num.samples, ncol = num.traits)
	for(i in 1:num.samples){
		traits.which <- sample(1:ncol(pheno), num.traits)
		sub.pheno <- pheno[,traits.which]
		sampled.traits[i,] <- colnames(pheno)[traits.which]
		total.ind[i] <- length(which(!is.na(rowSums(sub.pheno))))
		}
		
	full.mat <- cbind(sampled.traits, total.ind)
	ordered.combos <- t(apply(full.mat, 1, function(x) c(sort(x[1:num.traits]), x[(num.traits+1)])))
	
	u_combos <- unique(ordered.combos)
	sorted.comb <- u_combos[order(as.numeric(u_combos[,ncol(u_combos)]), decreasing = TRUE),]

	return(sorted.comb)
	
	
}