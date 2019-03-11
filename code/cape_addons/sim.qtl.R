#This function adds a QTL to a data object
#specifically for testing kinship corrections
#the added QTL is based on population structure
#allele frequency is the frequency of the causative
#allele in the different populations, not minor allele
#frequency. allele.freq must be the same length as the
#number of populations in the data.obj. Alleles will
#be simulated such that the different populations
#have the given frequency of the causative allele


sim.qtl <- function(data.obj, pop.assig, allele.freq, effect.size = 1, pheno = 1, qtl.position = NULL){
	
	
	
	#simulate genotypes for the different populations
	#in the data object
	num.pop <- length(unique(pop.assig))	
	if(is.null(dim(allele.freq))){
		allele.freq <- as.matrix(allele.freq, nrow = 1)
		}
	if(is.null(dim(effect.size))){
		effect.size <- as.matrix(effect.size, ncol = 1)
		}
	if(dim(allele.freq)[2] != num.pop){
		stop("allele.freq must be a vector with length equal to the number of populations")
		}
		
	allele.mat <- matrix(NA, nrow = length(pop.assig), ncol = dim(allele.freq)[1])
	for(p in 1:num.pop){ #for each distinct population in the data set
		pop.locale <- which(pop.assig == p)
		for(m in 1:dim(allele.freq)[1]){ #for each marker we are simulating
			allele.mat[pop.locale,m] <- sim.marker(length(pop.locale), allele.freq[m,p])
			}
		}
	#Now we alter the phenotype to be affected by the genotype
	new.pheno <- data.obj$pheno[,pheno]
	qtl.pheno <- new.pheno + allele.mat %*% effect.size


	model <- lm(qtl.pheno~allele.mat)
	print(summary(model))

		
	#pick random places to insert the markers
	if(is.null(qtl.position)){
		rnd.pos <- sort(sample(1:dim(data.obj$geno)[2], dim(allele.mat)[2], replace = FALSE))
		}else{
		rnd.pos <- qtl.position	
		}
	for(i in 1:length(rnd.pos)){
		data.obj$geno[,rnd.pos[i]] <- allele.mat[,i]
		}
	data.obj$pheno[,pheno] <- qtl.pheno
	
	orig.colnames <- colnames(data.obj$geno.for.pairscan)
	data.obj$geno.for.pairscan <- cbind(data.obj$geno.for.pairscan, data.obj$geno[,rnd.pos,drop=FALSE]) #add this for marking true QTL in the singlescan
	
	colnames(data.obj$geno.for.pairscan) <- c(orig.colnames, data.obj$marker.names[rnd.pos])
	
	
	return(data.obj)
	
}