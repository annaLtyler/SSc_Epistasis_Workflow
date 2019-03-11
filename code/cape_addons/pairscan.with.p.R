#This script performs the pairwise scan on all markers
#It takes in the data as a cross object.
#The user has the choice to scan the eigentraits (default)
#or the original phenotypes.
#This script also calls the function to do permutations
#on the 2D scan. It adds the genome-wide threshold for
#the 2D scan to the data object
#n.top.markers is used in generating the null. A permutation
#of the singlescan is run, and the n top markers are used
#in a permutation of the pairscan. if n.top.markers is NULL,
#it defaults to the number of markers in geno.for.pairscan

#scan.what = "eigentraits"; n.perm = 10; min.per.genotype = 6; use.pairs.threshold = TRUE; verbose = TRUE; num.pairs.limit = 1e4; num.perm.limit = 1e7


pairscan.with.p <- function(data.obj, kinship.obj, scan.what = c("eigentraits", "raw.traits"), n.perm = NULL, min.per.genotype = NULL, max.pair.cor = NULL, n.top.markers = NULL, verbose = FALSE, num.pairs.limit = 1e4, num.perm.limit = 1e7) {

	if(is.null(n.perm)){
		stop("The number of permutations must be specified.")
		}

	#If the user does not specify a scan.what, 
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentraits, otherwise, use raw phenotypes
	type.choice <- c(grep("eig", scan.what), grep("ET", scan.what), grep("et", scan.what)) #look for any version of eigen or eigentrait, the user might use.
	if(length(type.choice) > 0){ #if we find any, use the eigentrait matrix
		pheno <- data.obj$ET
		if(is.null(pheno)){stop("There are no eigentraits. Please set scan.what to raw.traits, or run get.eigentraits().")}
		}else{
			pheno <- data.obj$pheno #otherwise, use the raw phenotype matrix
			}

	num.pheno <- dim(pheno)[2]
	pheno.names <- colnames(pheno)
	geno <- data.obj$geno.for.pairscan
	covar.flags <- data.obj$covar.for.pairscan


	if(is.null(geno)){
		stop("select.markers.for.pairscan() must be run before pairscan()")
		}
		
		
	num.markers <- dim(geno)[2]
	
	#fill in a matrix to index the marker pairs
	marker.matrix <- pair.matrix(colnames(geno))
	pared.marker.mat <- get.pairs.for.pairscan(geno, min.per.genotype = min.per.genotype, max.pair.cor = max.pair.cor, verbose = verbose)

	num.pairs <- dim(pared.marker.mat)[1]
	
	if(num.pairs == 0){
		stop("There are no pairs to test. Try lowering min.per.genotype or raising max.pair.cor.")
		}

	if(!is.null(num.pairs.limit) && num.pairs > num.pairs.limit){
		cat("\nThe number of pairs (",num.pairs,") exceeds ", num.pairs.limit, ".\n", sep = "")
		go.on <- readline(prompt = "Do you want to continue (y/n)?\n")
		if(length(grep("n", go.on))){
			message("Stopping pairwise scan...\n")
			return(pairscan.obj)
		}else{
			cat("Continuing pairwise scan...\n")
		}
	}

	
	
	#make a list to hold the results. 
	#Tables from each of the phenotypes will be
	#put into this list
	results.list <- list()
	results.perm.list <- list()


	#run one.pairscan for each phenotype with results in scanone.result
	# if n.perm > 1 then call one.pairscan.perm()
	if(missing(kinship.obj)){
		for(p in 1:num.pheno){ 
	
			if(verbose){
				cat("\nScanning phenotype ", colnames(pheno)[p], ":\n", sep = "")
				}
			
			pairscan.results <- one.pairscan.with.p(phenotype.vector = pheno[,p], genotype.matrix = geno, covar.vector = covar.flags[,p], pairs.matrix = pared.marker.mat, n.perm = 0, verbose = verbose)
			results.list[[p]] <- pairscan.results[[1]]
			# results.perm.list[[p]] <- pairscan.results[[2]]
			} #end looping over phenotypes
		}else{
		#if we have a kinship object, we do one marker pair at a time
		#and pass in the corrected values

		chr.pairs <- strsplit(names(kinship.obj[[1]]), ",")
		chr.pair.mat <- cbind(unlist(lapply(chr.pairs, function(x) x[1])), unlist(lapply(chr.pairs, function(x) x[2])))
		
			for(p in 1:num.pheno){
				if(verbose){
				cat("\nScanning phenotype ", colnames(pheno)[p], ":\n", sep = "")
				}

				covar.vector <- covar.flags[,p]
				num.covar <- length(which(covar.vector == 1))
				
				effects.mat <- matrix(NA, nrow = dim(pared.marker.mat)[1], ncol = num.markers+num.covar+5)
				se.mat <- matrix(NA, nrow = dim(pared.marker.mat)[1], ncol = num.markers+num.covar+5)
				cov.mat <- matrix(NA, nrow = dim(pared.marker.mat)[1], ncol = 9)
				
				#pull out the list of corrections for the current phenotype
				kin <- kinship.obj[[p]]

				for(m in 1:dim(pared.marker.mat)[1]){
					if(verbose){report.progress(m, dim(pared.marker.mat)[1])}
					
					marker.locale <- match(pared.marker.mat[m,], colnames(data.obj$geno))
					marker.chr <- matrix(data.obj$chromosome[marker.locale], nrow = 1)
					chr.locale <- intersect(which(chr.pair.mat[,1] == marker.chr[1]), which(chr.pair.mat[,2] == marker.chr[2]))
								
					kindat <- kin[[chr.locale]]
					
					new.pheno <- kindat$corrected.pheno
					new.geno <- kindat$corrected.geno
					err.cov <- kindat$err.cov

					int.term = solve(err.cov) %*% new.geno[,marker.locale[1]]*new.geno[,marker.locale[2]]
		           
					pairscan.results <- one.pairscan(phenotype.vector = new.pheno, genotype.matrix = new.geno, int = int.term, covar.vector = covar.flags[,p], pairs.matrix = pared.marker.mat[m,,drop=FALSE], n.perm = 0, verbose = FALSE)
	
					effects.mat[m,] <- pairscan.results[[1]]$pairscan.effects
					se.mat[m,] <- pairscan.results[[1]]$pairscan.se
					cov.mat[m,] <- pairscan.results[[1]]$model.covariance
					} #end looping over marker pairs
					colnames(effects.mat) <- colnames(pairscan.results[[1]]$pairscan.effects)
					colnames(se.mat) <- colnames(pairscan.results[[1]]$pairscan.se)
										
					#add the results for the phenotype
					pheno.results <- list(effects.mat, se.mat, cov.mat)
					names(pheno.results) <- c("pairscan.effects", "pairscan.se", "model.covariance")
					results.list[[p]] <- pheno.results
				}	#end looping over phenotypes	
			}
	 
	 
	#generate the null distribution using pairscan.null
	num.pairs <- choose(dim(data.obj$geno.for.pairscan)[2], 2)
	total.perm = n.perm*num.pairs

	#calculate the number of top markers to use based
	#on the thresholding of the singlescan
	if(is.null(n.top.markers)){
		n.top.markers <- dim(data.obj$geno.for.pairscan)[2]
		}

	pairscan.obj <- pairscan.null(data.obj, scan.what = scan.what, total.perm = total.perm, n.top.markers = n.top.markers, verbose = verbose, min.per.genotype = min.per.genotype, max.pair.cor = max.pair.cor)


	names(results.list) <- pheno.names


	pairscan.obj$"pairscan.results" <- results.list	  #add the results to the data object



	return(pairscan.obj) #and return it
	
}
