#This script plots the histogram for var to pheno effects
#on a given chromosome.
#shade.over indicates where you would like to shade the 
#empirical distribution.
#You can shade it over a given percentile of the null 
#distribution. The default is to shade the portion of
#the empirical distribution over the 95th percentil of
#the null distribution.
#You can also shade based on a Holm's adjusted p value
#To do this, set shade.over to c.pval and set per.or.p
#to the adjusted p value. This will show any part of
#the distribution that is above the effect level corresponding
#to the adjusted p value. 


plot.effects.dist <- function(data.obj, perm.obj, chr = NULL, marker = NULL, shade.over = c("c.pval", "percentile"), perc.or.p = 0.05, pval.correction = c("holm", "fdr", "lfdr")){
	
	library(fdrtool)
	
	if(length(pval.correction) > 1){
		pval.correction = "fdr"
		}
	dir.inf <- perm.obj$var.to.pheno.test.stat
	dir.inf.perm <- perm.obj$var.to.pheno.test.stat.perm
	all.pheno <- names(dir.inf)
	
	if(is.null(dir.inf)){
		stop("direct.influence() must be run before plotting the effects")
		}

	if(is.null(chr) && is.null(marker)){
		stop("Either a marker or a chromosome must be specified.")
		}
	if(!is.null(chr) && !is.null(marker)){
		stop("Please specify either a marker or a chromosome. Not both.")
		}

	if(length(shade.over) == 1 && shade.over == "none"){
		shade = FALSE
		}else{
		shade = TRUE	
		}

	if(length(grep("pval", shade.over)) > 0){
		shade.above.percentile <- FALSE
		}else{
		shade.above.percentile <- TRUE
		}


	if(!is.null(chr)){
		#test to see if we are looking for a region
		region.test <- is.character(chr)
		
		if(region.test){
			region.locale <- c(which(names(data.obj $linkage.blocks.collapsed) == chr), which(names(data.obj $linkage.blocks.collapsed) == paste("Chr", chr, sep = "")))
			if(length(region.locale) == 0){
				stop("I can't find that chromosomal region. Please check the spelling and that linkage.blocks() has been run.")
				}
			marker.locale <- data.obj$linkage.blocks.collapsed[[region.locale]]
			}
		
		if(!region.test){
			chrom <- data.obj$chromosome	
			marker.locale <- which(chrom %in% chr)
			if(length(marker.locale) == 0){
				stop("I couldn't find any markers on chromosome", chr)
				}
			}

		marker.name <- marker
		if(length(dim(data.obj$geno)) == 2){
			marker <- colnames(data.obj$geno)[marker.locale]
			}else{
			marker <- data.obj$marker.names[marker.locale]
			}
		}else{
		marker.name <- marker
		if(length(dim(data.obj$geno)) == 2){
			marker <- colnames(data.obj$geno)[which(data.obj$marker.names %in% marker)]
			}else{
			locus.dim <- which(names(dimnames(data.obj$geno)) == "locus")
			marker <- dimnames(data.obj$geno)[[locus.dim]][which(data.obj$marker.names %in% marker)]	
			}
		}
		
		layout.mat <- get.layout.mat(length(all.pheno), "landscape")
		
		dev.new(width = dim(layout.mat)[2]*3, height = dim(layout.mat)[1]*4)
		layout(layout.mat)
		for(ph in 1:length(all.pheno)){
			pheno.dir.inf <- dir.inf[[ph]]
			pheno.dir.inf.perm <- dir.inf.perm[[ph]]

			if(length(dim(data.obj$geno)) == 2){
				all.dir.inf.locale <- which(pheno.dir.inf[,1] %in% marker)
				}else{
				just.markers <- sapply(strsplit(rownames(pheno.dir.inf), "_"), function(x) x[1])
				all.dir.inf.locale <- which(just.markers %in% marker)	
				}
			if(length(all.dir.inf.locale) == 0){
				dev.off()
				if(!is.null(chr)){
					stop("I can't find markers on chr ", chr,". Perhaps they weren't tested in the pair wise scan.", sep = "")				
					}else{
					stop("I can't find the marker ", marker,". Perhaps it wasn't tested in the pair wise scan.", sep = "")
					}
				}
			# all.dir.inf.perm.locale <- which(pheno.dir.inf.perm[,1] %in% marker)
			# all.dir.inf.perm.locale <- 1:dim(pheno.dir.inf.perm)[1]

			# null.dist <- abs(as.numeric(pheno.dir.inf.perm[,"t.stat"]))
			# emp.dist <- abs(as.numeric(pheno.dir.inf[all.dir.inf.locale,"t.stat"]))

			null.dist <- as.numeric(pheno.dir.inf.perm[,"t.stat"])
			emp.dist <- as.numeric(pheno.dir.inf[all.dir.inf.locale,"t.stat"])
			

			if(length(emp.dist) > 1){
				xmin <- min(c(density(null.dist)$x, density(emp.dist)$x)); xmax <- max(c(density(null.dist)$x, density(emp.dist)$x))
				ymax <- max(c(density(null.dist)$y, density(emp.dist)$y))*1.15
				}else{
				xmin <- min(c(null.dist, emp.dist)); xmax <- max(c(null.dist, emp.dist))
				ymax <- max(c(density(null.dist)$y))*1.15
				}

			if(!is.null(chr)){
				chr.text <- paste(chr, collapse = ", ")
				main.label <- paste("Chromosome", chr.text, "\nPhenotype:", all.pheno[ph])
				}else{
				marker.text <- paste(marker.name, collapse = ", ")
				main.label <- paste("Marker", marker.text, "\nPhenotype:", all.pheno[ph])
				}
			plot(density(null.dist), col = "red", xlim  = c(xmin, xmax), xlab = "Direct Influence t Statistic", main = main.label, ylim = c(0, ymax))
			if(length(emp.dist) > 1){
				points(density(emp.dist), col = "blue", type = "l")
				}else{
				arrows(emp.dist, ymax*0.1, emp.dist, 0, col = "blue", lwd = 3, length = 0.1)
				}

			#add a shaded region to the emp dist above the percentile in the null dist
		if(shade){
			if(shade.above.percentile){
				perc <- as.vector(quantile(null.dist, perc.or.p/100))
				}else{
				#get the p values of the distribution
				all.p <- apply(matrix(emp.dist, ncol = 1), 1, function(x) length(which(null.dist >= x))/length(null.dist))
				if(pval.correction == "holm"){
					c.p <- p.adjust(all.p, method = "holm")
					}
				if(pval.correction == "fdr"){
					fdr.out <- suppressWarnings(fdrtool(all.p, statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "pct0"))
					c.p <- fdr.out$qval
					}
				if(pval.correction == "lfdr"){
					fdr.out <- suppressWarnings(fdrtool(all.p, statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "pct0"))
					c.p <- fdr.out$lfdr
					}
					
				sig.locale <- which(c.p <= perc.or.p)
				if(length(sig.locale) > 0){
					perc <- emp.dist[min(sig.locale)]
					}else{
						perc <- max(null.dist)
						}
				}
			
			
			if(length(emp.dist) > 1){
				emp.dist.greater <- which(density(emp.dist)$x >= perc)
				}else{
				emp.dist.greater <- which(emp.dist >= perc)	
				}
				
			
						
			if(length(emp.dist.greater) > 0){
				x1 <- min(emp.dist.greater)
				x2 <- max(emp.dist.greater)
				if(length(emp.dist) > 1){
					polygon(x = density(emp.dist)$x[c(x1, x1:x2, x2)], y = c(0, density(emp.dist)$y[x1:x2], 0), col = "blue")
					}	
						
				}
			}

			legend("topright", legend = c("Null", "Observed"), lty = 1, col = c("red", "blue"), cex = 0.7)
			}
			
	
		
	
}