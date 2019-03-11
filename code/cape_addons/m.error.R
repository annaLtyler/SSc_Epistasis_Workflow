#This function looks at the properties of 
#markers with large standardized m values
#and very small std errors

m.error <- function(data.obj, pairscan.obj, geno.obj){
	require(timeDate)
	require(GGally)
	
	if(missing(pairscan.obj)){
		pairscan.results <- data.obj$pairscan.results	
		pairscan.perm <- data.obj$pairscan.perm		
		if(is.null(pairscan.results)){
			stop("I cannot find the pairscan results in the data object.")
			}
		}else{
		pairscan.results <- pairscan.obj$pairscan.results	
		pairscan.perm <- pairscan.obj$pairscan.perm
		}
		if(is.null(pairscan.results)){
			stop("I cannot find the pairscan results.")
			}		

	geno <- get.geno(data.obj, geno.obj)
	# geno[which(geno == 0.5)] <- 1
	#====================================================================================================
	# extract all the m values and their errors
	#====================================================================================================
	var.inf <- data.obj$var.to.var.influences
	var.inf.perm <- data.obj$var.to.var.influences.perm
	
	null.m12 <- as.numeric(var.inf.perm[,"m12"])
	null.m21 <- as.numeric(var.inf.perm[,"m21"])

	true.m12 <- as.numeric(var.inf[,"m12"])
	true.m21 <- as.numeric(var.inf[,"m21"])

	null.error.m12 <- as.numeric(var.inf.perm[,"m12.std.dev"])
	null.error.m21 <- as.numeric(var.inf.perm[,"m21.std.dev"])
	
	true.error.m12 <- as.numeric(var.inf[,"m12.std.dev"])
	true.error.m21 <- as.numeric(var.inf[,"m21.std.dev"])
	
	plot(log10(abs(null.m12))*sign(null.m12), null.m12/null.error.m12, xlab = "unstandardized", ylab = "standardized")
	plot(null.m21, null.m21/null.error.m21)

	#====================================================================================================
	# internal functions
	#====================================================================================================
	subdivide.markers <- function(effect.sizes, errors){
		std.effects <- effect.sizes/errors
		large.std <- which(abs(std.effects) > 1.5)
		small.std <- which(abs(std.effects) < 1.5)
		large.error <- which(errors > 1)
		small.error <- which(errors < 1)
		
		large.effect.large.error <- intersect(large.std, large.error)
		large.effect.small.error <- intersect(large.std, small.error)
		small.effect.large.error <- intersect(small.std, large.error)
		small.effect.small.error <- intersect(small.std, small.error)
		
		final.result <- list(large.effect.large.error, large.effect.small.error, small.effect.large.error, small.effect.small.error)
		names(final.result) <- c("large.std.effect_large.error", "large.std.effect_small.error", "small.std.effect_large.error", "small.std.effect_small.error")
		return(final.result)	
		}

	plot.group.effects <- function(sub.groups, perm = FALSE){
		if(perm){
			data.set <- var.inf.perm
			}else{
			data.set <- var.inf
			}

		par(mfrow = c(2,4))			
		all.effects <- lapply(sub.groups, function(x) log10(abs(as.numeric(data.set[x, "m12"])))*sign(as.numeric(data.set[x, "m12"])))
		xlim <- c(min(unlist(all.effects)), max(unlist(all.effects)))
		for(i in 1:length(sub.groups)){
			hist(all.effects[[i]], breaks = 100, main = names(sub.groups)[i], xlim = xlim, xlab = "log10 m12")
			}

		all.effects <- lapply(sub.groups, function(x) log10(abs(as.numeric(data.set[x, "m12"])))*sign(as.numeric(data.set[x, "m21"])))
		xlim <- c(min(unlist(all.effects)), max(unlist(all.effects)))
		for(i in 1:length(sub.groups)){
			hist(all.effects[[i]], breaks = 100, main = names(sub.groups)[i], xlim = xlim, xlab = "log10 m21")
			}
		mtext("Unstandardized log10 m12 and m21", line = -1.5, outer = TRUE)
		}
		
	get.maf <- function(marker.names){
		just.marker <- unlist(lapply(strsplit(marker.names, "_"), function(x) x[1]))
		marker.locale <- match(just.marker, dimnames(geno)[[3]])
		maf <- colMeans(geno[,2,marker.locale], na.rm = TRUE)
		return(maf)
		}
		
	plot.maf <- function(sub.groups, maf){
		par(mfrow = c(2,2))
		for(i in 1:length(sub.groups)){
			group.maf <- maf[sub.groups[[i]]]
			high.maf <- which(group.maf > 0.5)
			group.maf[high.maf] <- 1-group.maf[high.maf]
			hist(group.maf, breaks = 100, main = names(sub.groups)[i], xlab = names(sub.groups)[i])
			}
		}
		
	get.pair.pt.dist <- function(sub.groups, perm = FALSE){
		if(perm){
			data.set <- var.inf.perm
			}else{
			data.set <- var.inf
			}
		all.pt.rep <- vector(mode = "list", length = length(sub.groups))
		names(all.pt.rep) <- names(sub.groups)
		for(i in 1:length(sub.groups)){
			report.progress(i, length(sub.groups))
			m1 <- unlist(lapply(strsplit(data.set[sub.groups[[i]],1], "_"), function(x) x[1]))
			m2 <- unlist(lapply(strsplit(data.set[sub.groups[[i]],2], "_"), function(x) x[1]))
			m1.locale <- match(m1, dimnames(geno)[[3]])
			m2.locale <- match(m2, dimnames(geno)[[3]])
			# head(cbind(dimnames(geno)[[3]][m1.locale], m1))
			pt.rep <- apply(cbind(m1.locale, m2.locale), 1, function(x) as.vector(table(geno[,2,x[1]], geno[,2,x[2]])))
			if(class(pt.rep) == "matrix"){
				all.pt.rep[[i]] <- t(pt.rep)
				}else{
				all.pt.rep[[i]] <- Reduce("rbind", pt.rep)		
				}
			}
		return(all.pt.rep)		
		}
		
	plot.pt.rep <- function(pt.rep){
		par(mfrow = c(2,2))
		xlim <- c(min(unlist(pt.rep)), max(unlist(pt.rep)))
		for(i in 1:length(pt.rep)){
			hist(pt.rep[[i]], breaks = 100, main = names(pt.rep)[i], xlab = names(pt.rep)[i], xlim = xlim)
			}
		}
		
	get.pheno.dist <- function(sub.groups, perm = FALSE){
		if(perm){
			data.set <- var.inf.perm
			}else{
			data.set <- var.inf
			}
		
		all.pheno.rep <- vector(mode = "list", length = length(sub.groups))
		names(all.pheno.rep) <- names(sub.groups)
		for(i in 1:length(sub.groups)){
			report.progress(i, length(sub.groups))
			m1 <- unlist(lapply(strsplit(data.set[sub.groups[[i]],1], "_"), function(x) x[1]))
			m2 <- unlist(lapply(strsplit(data.set[sub.groups[[i]],2], "_"), function(x) x[1]))
			m1.locale <- match(m1, dimnames(geno)[[3]])
			m2.locale <- match(m2, dimnames(geno)[[3]])

			na.locale <- which(is.na(cbind(m1.locale, m2.locale)), arr.ind = TRUE)
			if(nrow(na.locale) > 0){
				m1.locale <- m1.locale[-na.locale[,1]]
				m2.locale <- m2.locale[-na.locale[,1]]
				}

			all.pheno.geno <- vector(mode = "list", length = length(m1.locale))
			
			for(g in 1:length(m1.locale)){
				geno.pairs <- cbind(geno[,2,m1.locale[g]], geno[,2,m2.locale[g]])
				geno.groups <- apply(genotypes, 1, function(x) intersect(which(geno.pairs[,1] == x[1]), which(geno.pairs[,2] == x[2])))
				pheno.mat <- matrix(NA, nrow = 4, ncol = ncol(data.obj$pheno))
				rownames(pheno.mat) <- apply(genotypes, 1, function(x) paste0(x[1], "_", x[2]))
				colnames(pheno.mat) <- colnames(data.obj$pheno)
				for(ph in 1:ncol(cross$pheno)){
					pheno.mat[,ph] <- unlist(lapply(geno.groups, function(x) mean(data.obj$pheno[x,ph])))
					}
				all.pheno.geno[[g]] <- pheno.mat
				}
			all.pheno.rep[[i]] <- all.pheno.geno
			}
		return(all.pheno.rep)
		}
		
	plot.pheno.rep <- function(pheno.rep){
		par(mfrow = c(2,2))
		for(i in 1:length(pheno.rep)){
			by.pheno <- Reduce("rbind", pheno.rep[[i]])
			all.dens <- apply(by.pheno, 2, density)
			xlim <- c(0,1)
			ylim <- c(min(unlist(lapply(all.dens, function(a) a$y))), max(unlist(lapply(all.dens, function(a) a$y))))
			cols <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
			plot.new()
			plot.window(xlim = xlim, ylim = ylim)
			for(j in 1:length(all.dens)){
				points(all.dens[[j]], type = "l", col = cols[j], lwd = 3)
				}	
			axis(1);axis(2)
			legend("topleft", col = cols[1:ncol(data.obj$pheno)], lty = 1, lwd = 3, legend = colnames(data.obj$pheno))
			mtext(names(pheno.rep)[i], side = 3, line = 1.5)
			}
		mtext("Phenotype Means By Group", side = 3, outer = TRUE, line = -1.5)
		}
	#====================================================================================================	

	#====================================================================================================
	# find the markers with large standardized effects
	# and further subdivide those into those with large and small errors
	#====================================================================================================
	m12.groups <- subdivide.markers(true.m12, true.error.m12)
	m21.groups <- subdivide.markers(true.m21, true.error.m21)
	m12.null.groups <- subdivide.markers(null.m12, null.error.m12)
	m21.null.groups <- subdivide.markers(null.m21, null.error.m21)

	#====================================================================================================
	# plot the unstandardized m values
	# for each of these groups of markers
	#====================================================================================================	
	plot.group.effects(m12.groups)
	plot.group.effects(m21.groups)	
	plot.group.effects(m12.null.groups)
	plot.group.effects(m21.null.groups)		

	#====================================================================================================	
	# look at the maf for these groups
	#====================================================================================================	
	m1.maf <- get.maf(var.inf[,1])
	m2.maf <- get.maf(var.inf[,2])
	m1.null.maf <- get.maf(var.inf.perm[,1])
	m2.null.maf <- get.maf(var.inf.perm[,2])
		
	plot.maf(m12.groups, m1.maf)
	plot.maf(m12.groups, m2.maf)	
	plot.maf(m21.groups, m1.maf)
	plot.maf(m21.groups, m2.maf)	

	plot.maf(m12.null.groups, m1.null.maf)
	plot.maf(m12.null.groups, m2.null.maf)	
	plot.maf(m21.null.groups, m1.null.maf)
	plot.maf(m21.null.groups, m2.null.maf)	
	
	#====================================================================================================	
	# look at the patient distribution by genotype for each group
	#====================================================================================================	
	true.pt.rep.by.m12 <- get.pair.pt.dist(sub.groups = m12.groups, perm = FALSE)
	true.pt.rep.by.m21 <- get.pair.pt.dist(m21.groups, FALSE)

	null.pt.rep.by.m12 <- get.pair.pt.dist(m12.null.groups, TRUE)
	null.pt.rep.by.m21 <- get.pair.pt.dist(m21.null.groups, TRUE)

	plot.pt.rep(true.pt.rep.by.m12)
	plot.pt.rep(true.pt.rep.by.m21)
	plot.pt.rep(null.pt.rep.by.m12)
	plot.pt.rep(null.pt.rep.by.m21)

	#====================================================================================================	
	# look at phenotype distribution of patient groups
	#====================================================================================================	
	genotypes <- pair.matrix(c(0,1), ordered = TRUE, self.pairs = TRUE)
	true.pheno.by.m12 <- get.pheno.dist(sub.groups = m12.groups, perm = FALSE)
	true.pheno.by.m21 <- get.pheno.dist(m21.groups, perm = FALSE)

	null.pheno.by.m12 <- get.pheno.dist(m12.null.groups, perm = FALSE)
	null.pheno.by.m21 <- get.pheno.dist(m21.null.groups, perm = FALSE)

	plot.pheno.rep(true.pheno.by.m12)
	plot.pheno.rep(true.pheno.by.m21)
	plot.pheno.rep(null.pheno.by.m12)
	plot.pheno.rep(null.pheno.by.m21)
	
	
	
}