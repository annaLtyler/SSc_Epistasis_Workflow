#This function plots the variance-covariance,
#and delta null distributions for a cross object


plot.covar.dist <- function(data.obj, perc.delta.to.plot = 99, verbose = TRUE){
	

	#====================================================================================
	#begin internal functions
	#====================================================================================
		get.beta.mat <- function(scan.two.results, marker.pair.number){
			#the beta matrix is composed of the coefficients from each pairwise
			#marker model (except for the interaction coefficient)
			#we use all markers in each row and set the non-covariate entries to 0
			num.pheno <- length(scan.two.results)
			beta.mat <- sapply(scan.two.results, function(x) as.numeric(x[[1]][marker.pair.number,(dim(scan.two.results[[1]][[1]])[2]-2):(dim(scan.two.results[[1]][[1]])[2])]))
			rownames(beta.mat) <- c("marker1", "marker2", "interaction")
			return(beta.mat)	
			}


		get.se.mat <- function(scan.two.results, marker.pair.number){
			#the beta matrix is composed of the coefficients from each pairwise
			#marker model (except for the interaction coefficient)
			#we use all markers in each row and set the non-covariate entries to 0
			num.pheno <- length(scan.two.results)
			se.mat <- sapply(scan.two.results, function(x) as.numeric(x[[2]][marker.pair.number,(dim(scan.two.results[[1]][[1]])[2]-2):(dim(scan.two.results[[1]][[1]])[2])]))
			rownames(se.mat) <- c("marker1", "marker2", "interaction")
			return(se.mat)	
			}


		
		get.cov.mat <- function(scan.two.results, marker.pair.number){
			#the variance-covariance matrix is block diagonal, and 
			#contains three times the number of rows and columns
			#as scanned traits. The blocks contain the variance-
			#covariance matrix from each trait for the marker pair
			cov.mat <- matrix(0, num.pheno*3, num.pheno*3)
			pheno.num <- 1:num.pheno
			start.row.col <- (pheno.num*3)-2
			end.row.col <- (pheno.num*3)
			for(ph in 1:num.pheno){
				cov.mat[start.row.col[ph]:end.row.col[ph], start.row.col[ph]:end.row.col[ph]] <- matrix(scan.two.results[[ph]][[3]][marker.pair.number,], ncol = 3)
				}	
			return(cov.mat)
			}
	#====================================================================================

		#check to see if phenotypes or eigentraits were scanned
		pheno.names <- names(data.obj$pairscan.results)	
		pheno.check <- match(pheno.names, colnames(data.obj$pheno))
		if(length(which(!is.na(pheno.check))) == 0){ #if we scanned eigentraits
			num.pheno <- dim(data.obj$ET)[2] #the number of phenotypes
			names.pheno <- colnames(data.obj$ET)
			}else{
			num.pheno <- dim(data.obj$pheno)[2] #the number of phenotypes
			names.pheno <- colnames(data.obj$pheno)
			}					


	### For all marker pairs calculate activity and IC
	n.gene <- dim(data.obj$geno.for.pairscan)[2] #the number of genes used in the pairscan

		if(is.null(data.obj$pairscan.perm)){
			stop("pairscan() with permutations must be run before error.prop()")
			}
		n.perm <- dim(data.obj$pairscan.perm[[1]][[1]])[1]/dim(data.obj$pairscan.results[[1]][[1]])[1]
		marker.mat <-  data.obj$pairscan.perm[[1]][[1]][,1:2]#get all the pairs that were tested in the pair scan
    	scan.two.results <- data.obj$pairscan.perm  #the coefficient matrices from the 2D scan			

   	colnames(marker.mat) <- c("marker1", "marker2")
	n.pairs <- length(marker.mat[,1]) #number of pairs of genes


	beta.main <- t(get.beta.mat(scan.two.results, 1)) ### Extract Main effect and interactions
	beta.se <- t(get.se.mat(scan.two.results, 1)) ### Extract Main effect and interactions
	beta.cov <- get.cov.mat(scan.two.results, 1) ### Extract Covars


	coeff.var <- matrix(NA, nrow = n.pairs, ncol = dim(beta.cov)[2])
	coeff.covar <- matrix(NA, nrow = n.pairs, ncol = length(as.vector(beta.cov[which(beta.cov != 0)])))
	delta.mat <- matrix(NA, nrow = n.pairs, ncol = 2)

	for(p in 1:length(marker.mat[,1])){

		if(verbose){
			report.progress(p, n.pairs, percent.text = 10, percent.dot = 2)
			}

		beta.main <- t(get.beta.mat(scan.two.results, p)) ### Extract Main effect and interactions
		beta.se <- t(get.se.mat(scan.two.results, p)) ### Extract Main effect and interactions
		beta.cov <- get.cov.mat(scan.two.results, p) ### Extract Covars
		non.zero <- which(beta.main[1:2,] != 0)
		if(length(non.zero) > 0){
			delta.mat[p,] <- get.delta(marker.mat[p,], beta.main, beta.se, beta.cov)
			}

		coeff.var[p,] <- diag(beta.cov)
		if(length(c(as.vector(beta.cov[which(beta.cov != 0)]))) == dim(coeff.covar)[2]){
			coeff.covar[p,] <- c(as.vector(beta.cov[which(beta.cov != 0)]))
			}
		}
		
	na.vals <- which(is.na(coeff.covar)) 
	if(length(na.vals) > 0){
		coeff.covar <- coeff.covar[-na.vals]
		}
		
	# all.results <- list(coeff.var, coeff.covar, delta.mat)
	# names(all.results) <- c("coeff.var", "coeff.covar", "delta.mat")
	# saveRDS(all.results, "delta.and.covar.matrices.RData")

	#find where the bulk of the deltas are
	test.seq <- seq(1, round(max(abs(delta.mat))), 10)
	prop.vals <- 0
	test.count <- 1
	while(prop.vals <= perc.delta.to.plot/100){
		prop.vals <- length(intersect(which(delta.mat > test.seq[test.count]*-1), which(delta.mat < test.seq[test.count])))/length(delta.mat)
		test.count <- test.count + 1
		# print(prop.vals)
		}
	x.lim <- c(test.seq[test.count+1]*-1, test.seq[test.count+1]) #go out to the next value to make sure we really capture the distribution
	
	pdf("Null.Covar.Matrix.pdf", width = 10, height = 4)
	par(mfrow = c(1,3))
	plot(density(coeff.var), main = "Variance Null", lwd = 3, cex.lab = 1.5, cex.axis = 2, cex.main = 2)
	plot(density(coeff.covar), main = "Covariance Null", lwd = 3, cex.lab = 1.5, cex.axis = 2, cex.main = 2)
	plot(density(delta.mat), main = "Delta Null", lwd = 3, cex.lab = 1.5, cex.axis = 2, cex.main = 2, xlim = x.lim)
	dev.off()
	
}