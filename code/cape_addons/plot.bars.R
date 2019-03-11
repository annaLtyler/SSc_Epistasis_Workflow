#This function plots the phenotypes from individual motifs
#for example purposes
#By default this function uses 0, 0.5, and 1 as genotype symbols,
#but you can specify alternate symbols with genotype.symbols.


plot.bars <- function(marker.geno, marker.name, phenotype.vals, pheno.name, ref.centered = FALSE, error.bars = "none", ymin = NULL, ymax = NULL, marker.pvals = marker.pvals, genotype.symbols = NULL, cex.axis = 1.5, text.cex = 1.5){

		error.bar.width = 0.1
		addline.width = 0.5
		addline.offset = 0.55
		error.bar.lwd = 2
		text.offset = 0.1
		
		if(dim(marker.geno)[2] == 1){
			genotypes <- sort(unique(marker.geno[,1]))
			pheno.vals <- lapply(genotypes, function(x) phenotype.vals[which(marker.geno[,1] == x)])
			pheno.means <- unlist(lapply(pheno.vals, function(x) mean(x, na.rm = TRUE)))
			if(error.bars == "sd"){
				pheno.error <- unlist(lapply(pheno.vals, function(x) sd(x, na.rm = TRUE)))
				}
			if(error.bars == "se"){
				pheno.error <- unlist(lapply(pheno.vals, function(x) sd(x, na.rm = TRUE)/sqrt(length(x))))
				}
			if(error.bars == "none"){
				pheno.error <- vector(mode = "list", length = length(genotypes))
				for(i in 1:length(pheno.error)){
					pheno.error[[i]] <- 0
					}
				}
			pheno.error[which(is.na(pheno.error))] <- 0
			if(is.null(ymin)){ymin <- min(c(pheno.means - pheno.error, 0))}
			if(is.null(ymax)){ymax <- max(pheno.means + pheno.error)}
			plot.height = ymax - ymin
				
			a <- barplot(pheno.means, ylim = c(ymin, ymax*1.1))
			segments(a, pheno.means-pheno.error, a, pheno.means+pheno.error, lwd = error.bar.lwd)
			abline(h = 0)
			#lower bar
			segments((a+error.bar.width), (pheno.means-pheno.error), (a-error.bar.width), pheno.means-pheno.error, lwd = error.bar.lwd)	
			#upper bar
			segments(a+error.bar.width, pheno.means+pheno.error, a-error.bar.width, pheno.means+pheno.error, lwd = error.bar.lwd)	
			if(is.null(genotype.symbols)){
				genotype.labels <- genotypes
				}else{
				genotype.labels <- genotype.symbols[1:length(genotypes)]
				}	
				
			text(x = genotypes, y = ymin-(plot.height*0.1), labels = genotype.labels)
			mtext(paste(pheno.name, "\np =",signif(marker.pvals, 2)), side = 3, line = 1.5)

			}else{ #instead if there are two markers

		marker1 <- marker.geno[,1]
		marker2 <- marker.geno[,2]
		
	#========================================================================
	# internal functions
	#========================================================================
	
	plot.bars.int <- function(pheno.vals, pheno.error, pheno.name, marker1, marker2, geno.n, glabels){
			
			if(length(pheno.vals) == 4){
				base.val = pheno.vals[1]
				}else{
				base.val = 0	
				}
			pred.add <- base.val + ((pheno.vals[2]-base.val) + (pheno.vals[3]-base.val))
			if(is.na(pred.add)){pred.add <- 0}
			pred.error <- pheno.error[2] + pheno.error[3]
	
			if(is.null(ymin)){
				if(error.bars != "none"){
					ymin <- min(c((pheno.vals-pheno.error), pred.add-pred.error, 0), na.rm = TRUE)
					}else{
					ymin <- min(c(pheno.vals, pred.add-pred.error, 0), na.rm = TRUE)	
					}
				}
			if(is.null(ymax)){
				if(error.bars != "none"){
					ymax <- max(c((pheno.vals+pheno.error), pred.add+pred.error), na.rm = TRUE)
					}else{
					ymax <- max(c(pheno.vals, pred.add+pred.error), na.rm = TRUE)	
					}
				}		
			
			
			full.range = ymax-ymin
			geno.n.y <- ymax+(full.range*0.1)
			ymax <- ymax*1.1
			
			label.y <- min(pheno.vals)-full.range*0.1
			a <- barplot(pheno.vals, ylim = c(ymin, ymax), axes = FALSE, xlim = c(0, length(pheno.vals)*1.5), names.arg = NA)
			# axis(2, at = signif(segment.region(ymin, ymax, 4, "ends"), 1), cex.axis = cex.axis) #this is an attempt to
					#get consecutive y axes to line up, but it didn't work
			axis(2, cex.axis = cex.axis)
			abline(h = 0)
			par(xpd = TRUE)

			text(x = 0, y = ymin-(full.range*0.1), labels = c(marker.name[1]), adj = 1, cex = text.cex)
			text(x = 0, y = ymin-(full.range*0.22), labels = c(marker.name[2]), adj = 1, cex = text.cex)
			text(x = c(a), y = ymin-(full.range*0.1), labels = c(glabels[1], glabels[2], glabels[1], glabels[2]), cex = text.cex)
			text(x = c(a), y = ymin-(full.range*0.22), labels = c(glabels[1], glabels[1], glabels[2], glabels[2]), cex = text.cex)
			if(error.bars != "none"){
				segments(a, pheno.vals-pheno.error, a, pheno.vals+pheno.error, lwd = error.bar.lwd)
				#lower bar
				segments((a+error.bar.width), (pheno.vals-pheno.error), (a-error.bar.width), pheno.vals-pheno.error, lwd = error.bar.lwd)	
				#upper bar
				segments(a+error.bar.width, pheno.vals+pheno.error, a-error.bar.width, pheno.vals+pheno.error, lwd = error.bar.lwd)	
				}
			segments(a[length(a)]+addline.width, pred.add, a[length(a)]-addline.width, pred.add, lty = 2, lwd = 2)
			poly.x <- c(a[length(a)]-addline.width, a[length(a)]+addline.width)
			poly.y <- c(pred.add-pred.error, pred.add+pred.error)
			polygon(x = c(poly.x, rev(poly.x)), y = rep(poly.y, each = 2), col = rgb(253/256,192/256,134/256, alpha = 0.5))
			arrows(a[length(a)]+addline.width+addline.offset, pred.add, a[length(a)]+addline.offset, pred.add, lty = 1, lwd = 2, length = 0.1)
			text(x = a[length(a)]+addline.width+addline.offset+text.offset, y = pred.add, labels = "Additive\nPrediction", adj = 0)
			text(x = a[,1], y = rep(geno.n.y, length(a[,1])), labels = geno.n)
			par(xpd = FALSE)
			mtext(pheno.name, side = 2, line = 2.5)
			}			

		#========================================================================

		marker1.geno <- marker1
		marker2.geno <- marker2
		
		#the genotype combinations used are based on the
		#maximum genotype present. We want to show the 
		#maximum effect for the "mutant" phenotype.
		min.geno <- min(c(marker1.geno, marker2.geno), na.rm = TRUE)
		max.geno <- max(c(marker1.geno, marker2.geno), na.rm = TRUE)
				
		if(is.null(genotype.symbols)){
			genotype.labels <- c(min.geno, max.geno)
			}else{
			genotype.labels <- genotype.symbols[1:2]	
			}
			w.w <- intersect(which(marker1.geno == min.geno), which(marker2.geno == min.geno))
			w.m <- intersect(which(marker1.geno == max.geno), which(marker2.geno == min.geno))
			m.w <- intersect(which(marker1.geno == min.geno), which(marker2.geno == max.geno))
			m.m <- intersect(which(marker1.geno == max.geno), which(marker2.geno == max.geno))
				
			pheno.list <- vector(mode = "list", length = 4)
			names(pheno.list) <- c(paste(min.geno, min.geno, sep = "/"), 
				paste(min.geno, max.geno, sep = "/"), 
				paste(max.geno, min.geno, sep = "/"), 
				paste(max.geno, max.geno, sep = "/"))
	
				pheno.list[[1]] <- phenotype.vals[w.w]
				pheno.list[[2]] <- phenotype.vals[w.m]
				pheno.list[[3]] <- phenotype.vals[m.w]
				pheno.list[[4]] <- phenotype.vals[m.m]
								
			
			geno.n <- unlist(lapply(pheno.list, length))
			pheno.means <- unlist(lapply(pheno.list, function(x) mean(x, na.rm = TRUE)))

			if(error.bars == "se"){
				pheno.error <- unlist(lapply(pheno.list, function(x) sd(x, na.rm = TRUE)/sqrt(length(x))))
				}
			if(error.bars == "sd"){
				pheno.error <- unlist(lapply(pheno.list, function(x) sd(x, na.rm = TRUE)))
				}

			if(!ref.centered){
				plot.bars.int(pheno.vals = pheno.means, pheno.error = pheno.error, pheno.name = pheno.name, marker1 = marker1, marker2 = marker2, geno.n = geno.n, glabels = genotype.labels)
				}else{
				#also make a plot using B6 as the reference point
				#center on B6 genotype
				ref.pheno <- pheno.means-pheno.means[1]
			
				pheno.means[which(is.na(pheno.means))] <- 0
				pheno.error[which(is.na(pheno.error))] <- 0
				plot.bars.int(pheno.vals = ref.pheno, pheno.error, pheno.name = pheno.name, marker1 = marker1, marker2 = marker2, geno.n = geno.n, glabel = genotype.labels)
				}		
			}
		} #end case for two markers