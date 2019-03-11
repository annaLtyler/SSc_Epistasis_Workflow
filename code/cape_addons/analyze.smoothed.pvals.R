#This function compares p values from smoothed few permutations
#to unsmoothed p values from large numbers of permutations
#num.perms <- seq(1000, 10000, 1000)

analyze.smoothed.pvals <- function(data.obj, num.perms = seq(100000, 1000000, 100000), num.resamples = 10){
		
	
	influences.org <- data.obj$var.to.var.influences
	influences.perm <- data.obj$var.to.var.influences.perm    
    	
    marker.mat <- influences.org[,1:2] #a matrix listing the names of all marker combinations
    colnames(marker.mat) <- c("marker1", "marker2")

	#### Combinine across permutations#####
	#get the t statistics for all permutations
	mat12.perm <- as.numeric(influences.perm[,3]) / as.numeric(influences.perm[,4])
	mat21.perm <- as.numeric(influences.perm[,5]) / as.numeric(influences.perm[,6])
	mat12.mat21.perm <- c(mat12.perm, mat21.perm)

	mat12 <- as.numeric(influences.org[,3]) / as.numeric(influences.org[,4])
	mat21 <- as.numeric(influences.org[,5]) / as.numeric(influences.org[,6])
	mat12.mat21 <- c(mat12, mat21)

	#find where the negative an positive mat12.mat21 values are
	neg.vals <- which(mat12.mat21 < median(mat12.mat21.perm))
	pos.vals <- which(mat12.mat21 > median(mat12.mat21.perm))
	colors <- rep(NA, length(mat12.mat21))
	colors[neg.vals] <- "blue"
	colors[pos.vals] <- "red"

	#========================================================================================================
	#internal functions
	#========================================================================================================
	
	#use above and below median
	get.emp.p <- function(emp.val, null.dist){
		if(emp.val < median(null.dist)){
			emp.p <- (length(which(null.dist <= emp.val))+1)/(length(which(null.dist <= median(null.dist)))+1)
			}else{
			emp.p <- (length(which(null.dist >= emp.val))+1)/(length(which(null.dist >= median(null.dist)))+1)
			}
		return(emp.p)
		}

		
	#two-sided p values
	two.sided.p <- function(emp.val, null.dist){
		emp.p <- length(which(abs(null.dist) >= abs(emp.val)))/length(null.dist)
		return(emp.p)
		}

	#========================================================================================================
	
	unsmoothed.p.mats <- vector(mode = "list", length = length(num.perms))
	smoothed.p.mats <- vector(mode = "list", length = length(num.perms))
	
	stage.count <- 1
	for(s in num.perms){
		#resample the adjusted null a given number of times and recalculate
		#pvalues using the standard 2-sided method
		cat("\nWorking on",s, "permutations.\n")
		smoothed.p.mat <- matrix(NA, nrow = length(mat12.mat21), ncol = num.resamples)
		unsmoothed.p.mat <- matrix(NA, nrow = length(mat12.mat21), ncol = num.resamples)
		
		pdf(paste("Smoothed.p.vs.asymmetric.p.", s, ".permutations.pdf", sep = ""))
		for(i in 1:num.resamples){
			report.progress(i, num.resamples)
			chopped.null <- sample(mat12.mat21.perm, s)

			#get p values in the skewed distribution using a one-sided calculation
			unsmoothed.pvals <- t(apply(matrix(mat12.mat21, ncol = 1), 1, function(x) get.emp.p(x, chopped.null)))
			unsmoothed.p.mat[,i] <- unsmoothed.pvals
		
			#get the p values for each marker pair from the same 
			#numbers of permutations with smoothed distributions
			m12.m21.adj <- symmetrize.null(null.dist = chopped.null, emp.dist = mat12.mat21, plot.dist = TRUE, "test.pdf")		
			adj.null <- m12.m21.adj$normalized.null
			adj.m12.m21 <- m12.m21.adj$normalized.empirical
			pvals <- apply(matrix(adj.m12.m21, ncol = 1), 1, function(x) two.sided.p(x, adj.null))
			smoothed.p.mat[,i] <- pvals
			}

			unsmoothed.p.mats[[stage.count]] <- unsmoothed.p.mat
			smoothed.p.mats[[stage.count]] <- smoothed.p.mat

			log.unsmoothed.p <- log10(unsmoothed.p.mat)*-1
			log.smoothed.p <- log10(smoothed.p.mat)*-1	
						
			smoothScatter(as.vector(log.unsmoothed.p), as.vector(log.smoothed.p), colramp = colorRampPalette(c("white", "black"),space = "rgb"), pch = 16, cex = 0.3, ylab = "-log10 p value from smoothed distribution", xlab = "one-sided -log10 p value from unsmoothed distribution", main = paste(s, "permutations"))
			abline(0,1, lwd = 2)
			dev.off()
			stage.count <- stage.count + 1
			}
		
		saveRDS(unsmoothed.p.mats, "unsmoothed.pmats.RData")
		saveRDS(smoothed.p.mats, "smoothed.pmats.RData")	

}