#plot interactions for given marker pair

plot.interaction <- function(data.obj, geno.obj = NULL, marker1, marker2, scan.what = c("ET", "normalized.trait", "raw.trait"), phenotype.name, covar = NULL, ymin = NULL, ymax = NULL, error.bars = "se", marker1.labels = NULL, marker2.labels = NULL){
	
	geno <- get.geno(data.obj, geno.obj)
	marker.vals <- get.marker.val(data.obj, geno.obj, c(marker1, marker2), alleles = "B")

	marker.geno <- apply(marker.vals, 2, bin.vector)

	genotypes <- sort(unique(c(marker.geno[,1], marker.geno[,2])))
	num.genotypes <- length(genotypes)
	
	
	pheno <- get.pheno(cross, scan.what = scan.what, covar)
	pheno.locale <- which(colnames(pheno) == phenotype.name)
	pheno.vals <- data.obj$pheno[,pheno.locale]

	errors <- get.interaction.error(marker.geno[,1], marker.geno[,2], pheno.vals, error.type = error.bars)	
			
					
			if(is.null(ymin)){
				ymin <- min(errors[[1]]-errors[[2]], na.rm = TRUE)
				}
			if(is.null(ymax)){
				ymax <- max(errors[[1]]+errors[[2]], na.rm = TRUE)
				}

			ylim <- c(ymin, ymax)
			
			not.na <- which(!is.na(pheno.vals))
			
			interaction.plot(marker.geno[not.na, 1], marker.geno[not.na, 2], pheno.vals[not.na], xlab = marker1, trace.label = marker2, ylab = "", lwd = 3, cex.lab = 2.5, main = phenotype.name, axes = FALSE, ylim = ylim, fixed = TRUE, cex.main = 2.5)

			plot.lim <- par("usr")
			min.area = plot.lim[1]*1.17; max.area = plot.lim[2]*0.85
			axis(1, at = segment.region(min.area, max.area, num.genotypes, "ends"), labels = genotypes, cex.axis = 2); axis(2, cex.axis = 2)
			
			if(error.bars != "none"){
				for(i in 1:length(errors$means[,1])){
					segments(1:length(colnames(errors$means)), (errors$means[i,]+errors[[2]][i,]), 1:length(colnames(errors$means)), (errors$means[i,]-errors[[2]][i,]))
					}
				}
		
	}
		