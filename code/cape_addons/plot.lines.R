#if 1 marker, plot mean phenotype values with error bars
#if 2 markers, plot interaction plot
# calculate error bars inside function
#marker.geno is a matrix of 1 or 2 columns
#marker.name is a vector of length 1 or 2
plot.lines <- function(data.obj, marker.geno, marker.name, phenotype.vals, phenotype.name, ymin = NULL, ymax = NULL, p.value, error.bars = "none", marker.pvals, ref.centered = FALSE){
	
	just.markers <- sapply(strsplit(marker.name, "_"), function(x) x[1])
	marker.chr <- rev(get.marker.chr(data.obj, just.markers))
	marker.chr[which(is.na(marker.chr))] <- 0

	chr.text <- unlist(lapply(marker.chr, function(x) if(x == 0){""}else{paste0("Chr", x, "_")}))

	error.bar.width = 0.15


	if(dim(marker.geno)[2] == 1){
		
		#bin the phenotype into the genotype values, which have already been binned
		geno.bins <- sort(unique(marker.geno[which(!is.na(marker.geno))]))

		pheno.bins <- vector(mode = "list", length = length(geno.bins))
		names(pheno.bins) <- geno.bins
			
		for(g in 1:length(geno.bins)){
			pheno.bins[[g]] <- phenotype.vals[which(marker.geno == geno.bins[g])]
			}

		#calculate mean and standard error for each group
		pheno.means <- unlist(lapply(pheno.bins, function(x) mean(x, na.rm = TRUE)))
		
		if(ref.centered){
			ref.val <- pheno.means[[1]]
			pheno.means <- pheno.means-ref.val
			pheno.bins <- lapply(pheno.bins, function(x) x - ref.val)
			}
		
		if(error.bars == "sd"){
			pheno.error <- unlist(lapply(pheno.bins, function(x) sd(x, na.rm = TRUE)))
			}
		if(error.bars == "se"){
			pheno.error <- unlist(lapply(pheno.bins, function(x) sd(x, na.rm = TRUE)/sqrt(length(x))))
			}
		if(error.bars == "none"){
			pheno.error <- rep(0, length(pheno.bins))
			}

		if(is.null(ymin)){
			ymin <- min(pheno.means-pheno.error, na.rm = TRUE)
			}
		if(is.null(ymax)){
			ymax <- max(pheno.means+pheno.error, na.rm = TRUE)
			}


		plot(pheno.means, type = "b", ylim = c(ymin, ymax), xlim = c(1,length(geno.bins)), axes = FALSE, ylab = "", xlab = paste0(chr.text, marker.name, sep = ""), main = paste(phenotype.name, "\np =", signif(marker.pvals, 2)), lwd = 3, cex.lab = 1.5, cex.main = 1.5)
		axis(1, at = 1:length(pheno.means), labels = geno.bins, cex.axis = 2)
		axis(2, cex.axis = 2)
		
		if(error.bars != "none"){
			segments(x0 = 1:length(pheno.means), y0 = pheno.means-pheno.error, y1 = pheno.means+pheno.error, lwd = 2)
			}
		}
		
		if(dim(marker.geno)[2] == 2){
			
			marker.geno.bins <- vector(mode = "list", length = 2)
			for(g in 1:length(marker.geno.bins)){
				marker.geno.bins[[g]] <- sort(unique(marker.geno[,g]))
				}
			
			errors <- get.interaction.error(marker.geno[,1], marker.geno[,2], phenotype.vals, error.type = error.bars)	
			
			if(ref.centered){
				ref.val <- errors$means[1,1]
				phenotype.vals <- phenotype.vals - ref.val
				errors$means <- errors$means - ref.val
				}	
					
			if(is.null(ymin)){
				all.errors <- errors[[2]]
				all.errors[which(is.na(all.errors))] <- 0
				ymin <- min(errors[[1]]-all.errors, na.rm = TRUE)
				}
			if(is.null(ymax)){
				all.errors <- errors[[2]]
				all.errors[which(is.na(all.errors))] <- 0
				ymax <- max(errors[[1]]+all.errors, na.rm = TRUE)
				}

			ylim <- c(ymin, ymax)
			
			not.na <- which(!is.na(phenotype.vals))
			
			interaction.plot(marker.geno[not.na,1], marker.geno[not.na,2], phenotype.vals[not.na], xlab = paste0(chr.text[2], marker.name[1]), trace.label = paste(chr.text[1], marker.name[2]), ylab = "", lwd = 3, cex.lab = 1.5, main = phenotype.name, axes = FALSE, ylim = ylim, fixed = TRUE, cex.main = 1.5)

			plot.lim <- par("usr")
			min.area = plot.lim[1]*1.17; max.area = plot.lim[2]*0.85
			axis(1, at = c(seq(min.area, max.area, (max.area-min.area)/(length(marker.geno.bins[[2]])-1))), labels = marker.geno.bins[[2]], cex.axis = 2); axis(2, cex.axis = 2)
			
			if(error.bars != "none"){
				for(i in 1:length(errors$means[,1])){
					segments(1:length(colnames(errors$means)), (errors$means[i,]+errors[[2]][i,]), 1:length(colnames(errors$means)), (errors$means[i,]-errors[[2]][i,]))
					}
				}
		}#end case for if there are two markers
		
	}#end function
		