#This function plots p values for m12/m21's with increasing numbers of
#permutations


plot.p.curve <- function(data.obj, min.perm = 1000, max.perm = "max", num.bins = 100){


	cur.dir <- getwd()
	r.dir <- "~/Documents/git_repositories/cape"
	setwd(r.dir)
	all.fun <- list.files(pattern = "*.R")
	for(i in 1:length(all.fun)){source(all.fun[i])}
	setwd(cur.dir)
	
	
	influences.org <- data.obj$var.to.var.influences
	influences.perm <- data.obj$var.to.var.influences.perm

	n.gene <- dim(data.obj$geno.for.pairscan)[2] #get the number of genes used in the pair scan
	n.pairs <- dim(data.obj$pairscan.results[[1]][[1]])[1] #the number of pairs scanned in the pairscan
    
    	
    marker.mat <- influences.org[,1:2] #a matrix listing the names of all marker combinations
    colnames(marker.mat) <- c("marker1", "marker2")

	#====================================================================================
	#internal functions
	#====================================================================================
	two.sided.p <- function(emp.val, null.dist){
		# emp.p <- length(which(abs(null.dist) >= abs(emp.val)))/length(null.dist)
		emp.p <- (length(which(abs(null.dist) >= abs(emp.val)))+1)/(length(null.dist)+1)
		return(emp.p)
		}
	
	one.sided.p <- function(emp.val, null.dist){
		if(emp.val < median(null.dist)){
			emp.p <- length(which(null.dist <= emp.val))/length(which(null.dist <= median(null.dist)))
			}else{
			emp.p <- length(which(null.dist >= emp.val))/length(which(null.dist >= median(null.dist)))
			}
		return(emp.p)
		}
	#====================================================================================
	
	

	#### Combinine across permutations#####
	#get the t statistics for all permutations
	mat12.perm <- as.numeric(influences.perm[,3]) / as.numeric(influences.perm[,4])
	mat21.perm <- as.numeric(influences.perm[,5]) / as.numeric(influences.perm[,6])
	mat12.mat21.perm <- c(mat12.perm, mat21.perm)

	if(max.perm == "max"){
		max.perm <- length(mat12.mat21.perm)
		}
	perm.seq <- round(segment.region(min.perm, max.perm, num.bins, "ends"))
	
	cat("Looking at", min.perm, "to", max.perm, "permutations.\n")	
	
	mat12 <- as.numeric(influences.org[,3]) / as.numeric(influences.org[,4])
	mat21 <- as.numeric(influences.org[,5]) / as.numeric(influences.org[,6])
	mat12.mat21 <- c(mat12, mat21)

	#look at the one-sided p values for skewed distributions
	cat("\nCalculating one-sided p values...\n")
	p.mat <- matrix(NA, nrow = length(mat12.mat21), ncol = length(perm.seq))
	for(i in 1:length(perm.seq)){
		report.progress(i, length(perm.seq))
		chopped.null <- mat12.mat21.perm[1:perm.seq[i]]
		pvals <- apply(matrix(mat12.mat21, ncol = 1), 1, function(x) one.sided.p(x, chopped.null))
		p.mat[,i] <- pvals
		}

	saveRDS(p.mat, "pval.seq.one.sided.RData")

	to.plot <- log10(p.mat)*-1
	inf.vals <- which(to.plot == Inf)
	if(length(inf.vals) > 0){
		to.plot[inf.vals] <- 10
		}

	pdf("pval.curves.one.sided.pdf")	
	ymin = min(to.plot); ymax = max(to.plot)
	plot(x = perm.seq, to.plot[1,], type = "l", ylim = c(ymin, ymax), ylab = "-log10 p value", main = paste("Unadjusted m12/m21 p values for individual marker pairs\nfrom", min(perm.seq), "to", max(perm.seq), "permutations"), xlab = "Permutations")
	for(i in 2:dim(p.mat)[1]){
		points(x = perm.seq, to.plot[i,], col = i, type = "l")
		}
	axis(1)
	axis(2)
	dev.off()



	#fix the skew empirically and look at the two-sided p value
	mat12.adj <- symmetrize.null(null.dist <- mat12.mat21.perm, emp.dist = mat12, plot.dist = FALSE, plot.title = "Transformation.m12.Null.pdf")
	mat21.adj <- symmetrize.null(null.dist <- mat12.mat21.perm, emp.dist = mat21, plot.dist = FALSE, plot.title = "Transformation.m21.Null.pdf")	

	adj.null <- mat12.adj$normalized.null
	adj.m12 <- mat12.adj$normalized.empirical
	adj.m21 <- mat21.adj$normalized.empirical
	adj.m12.m21 <- c(adj.m12, adj.m21)
 
 	cat("\nCalculating two-sided p values...\n")
	p.mat <- matrix(NA, nrow = length(adj.m12.m21), ncol = length(perm.seq))
	for(i in 1:length(perm.seq)){
		report.progress(i, length(perm.seq))
		chopped.null <- adj.null[1:perm.seq[i]]
		pvals <- apply(matrix(adj.m12.m21, ncol = 1), 1, function(x) two.sided.p(x, chopped.null))
		p.mat[,i] <- pvals
		}


	saveRDS(p.mat, "pval.seq.two.sided.RData")
	
	to.plot <- log10(p.mat)*-1
	inf.vals <- which(to.plot == Inf)
	if(length(inf.vals) > 0){
		to.plot[inf.vals] <- 10
		}

	pdf("pval.curves.two.sided.pdf")	
	ymin = min(to.plot); ymax = max(to.plot)
	plot(x = perm.seq, to.plot[1,], type = "l", ylim = c(ymin, ymax), ylab = "-log10 p value", main = paste("Unadjusted m12/m21 p values for individual marker pairs\nfrom", min(perm.seq), "to", max(perm.seq), "permutations"), xlab = "Permutations")
	for(i in 2:dim(p.mat)[1]){
		points(x = perm.seq, to.plot[i,], col = i, type = "l")
		}
	axis(1)
	axis(2)
	dev.off()

	
	
}