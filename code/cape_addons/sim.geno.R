sim.geno <- function(n.ind = 100, n.markers = 100, num.pop = 4, noise.sd = 0.1){
	

	# genotypes <- c(-1, 0, 1)
	genotypes <- c(0, 0.5, 1)

	set.genotypes2 <- function(marker.freq, pop.ind){
		p <- 1-marker.freq
		q <- marker.freq
		maj.hom <- p^2
		het <- 2*p*q
		min.hom <- q^2

		samp.genos <- sample(genotypes, pop.ind, prob = c(maj.hom, het, min.hom), replace = TRUE)
			
		return(samp.genos)
		}


	#generate an allele frequency for the overall 
	#population. Adjust this slightly to make separate
	#populations.
	allele.freq <- runif(n.markers, 0, 0.5)
		
	indV <- 1:n.ind
	pop.ind <- vector(mode = "list", length = num.pop)
	geno.mat <- matrix(NA, nrow = n.ind, ncol = n.markers)
	for(p in 1:num.pop){
		if(p == num.pop){
			pop.sample <- 1:length(indV)
			}else{
			pop.sample <- sort(sample(1:length(indV), round(n.ind/num.pop)))	
			}	
		pop.ind[[p]] <- indV[pop.sample]
		#adjust the allele frequencies for the population
		pop.allele.freq <- abs(allele.freq + rnorm(n.markers, 0, noise.sd))
		pop.genos <- lapply(pop.allele.freq, function(x) set.genotypes2(x, length(pop.sample)))
		pop.mat <- matrix(unlist(pop.genos), ncol = n.markers, byrow = FALSE)
		geno.mat[indV[pop.sample],] <- pop.mat
		indV <- indV[-pop.sample]
		}
	
	pop.assig <- rep(NA, n.ind)
	for(i in 1:length(pop.ind)){
		pop.assig[pop.ind[[i]]] <- i
		}
		
	rownames(geno.mat) <- 1:dim(geno.mat)[1]
	colnames(geno.mat) <- 1:dim(geno.mat)[2]
	
	result <- list("geno" = geno.mat, "pop" = pop.assig)	
	return(result)
	
		
}