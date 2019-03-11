#This script plots the final main effects 
#in the same style as the 1D scan


plotFinalMainEffects <- function(data.obj, chr = NULL, traits = NULL, standardized = TRUE, mark.chr = TRUE, plot.type = "l", overlay = FALSE, trait.colors = NULL, show.peak.markers = FALSE, mark.pvals = c(0.01, 0.05), show.linkage.blocks = FALSE, chr.label.cex = 0.5, block.delim.lwd = 1, block.delim.col = "blue", block.delim.lty = 3, pheno.labels = NULL){
	
	mark.block.boundaries <- FALSE
	D1.results <- data.obj$max.var.to.pheno.influence
	#order the results by marker order
	D1.results <- lapply(D1.results, function(x) x[order(x[,1]),])
	marker.names <- data.obj$marker.names
	ind.markers <- data.obj$geno.for.pairscan
		
	if(is.null(D1.results)){
		stop("direct.influence() must be run before plotting the results")
		}


	if(is.null(chr)){
		chr <- unique(data.obj$chromosome[which(colnames(data.obj$geno) %in% D1.results[[1]][,1])])
		}

	if(is.null(traits)){
		traits <- names(D1.results)
		}
		

	if(show.linkage.blocks){
		blocks <- data.obj$linkage.blocks.collapsed
		if(is.null(blocks)){
			stop("linkage.blocks must be run before linkage blocks can be plotted.")
			}
		if(!is.null(chr)){
			block.locale <- which(sapply(strsplit(names(blocks), "_"), function(x) x[1]) %in% paste("Chr", chr, sep = ""))
			just.blocks <- vector(mode = "list", length = length(block.locale))
			names(just.blocks) <- names(blocks)[block.locale]
			for(i in 1:length(just.blocks)){
				just.blocks[[i]] <- blocks[[block.locale[i]]]
				}
			blocks = just.blocks
			}
		}
	
	if(show.peak.markers){
		peaks <- data.obj$peak.markers
		}

	covar.flags <- data.obj$covar.flags
	col.mat <- matrix(NA, nrow = dim(covar.flags)[1], ncol = length(D1.results))

	if(!overlay){
		col.mat[covar.flags == 0] <- "black"
		col.mat[covar.flags == 1] <- "black"
		}else{
		if(is.null(trait.colors)){
			trait.colors <- c("black", "blue", "purple", "darkgreen")
			}
		if(length(trait.colors) < length(traits)){
		 	trait.colors <- rep(trait.colors, length(traits)/4)
		 	trait.colors <- trait.colors[1:length(traits)]
			}
		for(i in 1:length(traits)){
			col.mat[,i] <- trait.colors[i]
			}
		}


	chr.locale <- which(data.obj$chromosome %in% chr)
	markers.used <- which(colnames(data.obj$geno) %in% D1.results[[1]][,1])
	markers.which <- intersect(chr.locale, markers.used)
	results.rows <- which(D1.results[[1]][,1] %in% colnames(data.obj$geno)[markers.which])
	results.el <- which(names(D1.results) %in% traits)
	results.to.plot <- NULL
	for(r in results.el){
		if(standardized){
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][results.rows,"|t.stat|"])
			}else{
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][results.rows,"coef"])
			}	
		}


	rownames(results.to.plot) <- D1.results[[1]][results.rows,1]
	final.cols <- col.mat[results.rows, results.el]
	colnames(results.to.plot) <- names(D1.results)[results.el]


	if(is.null(pheno.labels)){
		pheno.labels <- colnames(results.to.plot)
		}
	
	#also pull out the corresponding peak markers
	if(show.peak.markers){
		peaks.to.plot <- vector(mode = "list", length = length(traits))
		names(peaks.to.plot) <- traits
		for(tr in 1:length(traits)){
			pheno.locale <- which(peaks[,"phenotype"] %in% traits[tr])
			trait.peaks <- peaks[pheno.locale,]
			chr.locale <- which(trait.peaks[,"marker.chr"] %in% chr)
			peaks.to.plot[[tr]] <- trait.peaks[chr.locale,]
			}
		}
	

	#find the effect sizes where the specified corrected p values lie
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(D1.results[[1]]))))
		sig.effect.levels <- NULL
		for(pv in 1:length(mark.pvals)){
			if(standardized){
				sig.effect.levels <- rbind(sig.effect.levels, unlist(lapply(D1.results, function(x) min(x[which(x[,var.sig.col] <= mark.pvals[pv]),"|t.stat|"]))))
				}else{
				sig.effect.levels <- rbind(sig.effect.levels, unlist(lapply(D1.results, function(x) min(x[which(x[,var.sig.col] <= mark.pvals[pv]),"coef"]))))	
				}
			}

	if(!overlay){
		layout.mat <- matrix(c(1:2), ncol = 1)
		}else{
		layout.mat <- matrix(1, 1, 1)
		}

	layout(layout.mat)	
	for(p in 1:length(results.to.plot[1,])){
		# dev.new(width = plot.width, height = plot.height)
		pheno.res <- results.to.plot[,p]
		if(show.peak.markers){
			peak.height <- rep(NA, length(pheno.res))
			peak.locale <- which(data.obj$marker.names %in% peaks.to.plot[[p]][,1])
			peak.height[peak.locale] <- pheno.res[peak.locale]*1.02
			}
		
		if(standardized){
			if(!overlay){
				all.vals <- c(pheno.res, 0)		
				}else{
				all.vals <- c(results.to.plot, 0)	
					}
			}else{
				if(!overlay){
					all.vals <- pheno.res
					}else{
					all.vals <- results.to.plot	
					}
			}

		#create the window
		if(p == 1 || !overlay){
			par(mar = c(3, 4, 3, 2) + 0.1)
			plot.new()
			plot.window(xlim = c(0, length(pheno.res)), ylim = c(min(all.vals), max(all.vals[is.finite(all.vals)])))
			
		
			#shade the chromosome regions
			if(mark.chr){
				markers.used.locale <- which(colnames(data.obj$geno) %in% rownames(results.to.plot))
				chr.id <- data.obj$chromosome[markers.used.locale]
				par(xpd = TRUE)
				for(ch in 1:length(chr)){
					x.min <- min(which(chr.id == chr[ch])); x.max <- max(which(chr.id == chr[ch]))
					if(ch %% 2 == 1){
						polygon(x = c(x.min, x.min, x.max, x.max), y = c(min(all.vals), max(all.vals), max(all.vals), min(all.vals)), col = "lightgray", border = NA)
						}
					if(chr[ch] == 0){
						text(x = x.max, y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = "Cov.", adj = 0, cex = chr.label.cex)
						}else{
						text(x = mean(c(x.min, x.max)), y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = chr[ch], cex = chr.label.cex)
						}
					}
				par(xpd = FALSE)
				}
		
				# axis(1, at = 1:length(pheno.res), labels = FALSE)
				abline(h = 0)
			    # lbl <- marker.names



			if(standardized){

				#add the marker labels to the x axis	    
			    # text(1:length(pheno.res), par("usr")[3] - 0.25, srt = 90, adj = 1, labels = NULL, xpd = TRUE, cex = 0.5)
				#add lines to indicate the significance thresholds
				for(pv in 1:length(mark.pvals)){
					abline(h = sig.effect.levels[pv,p], lty = pv)
					par(xpd = TRUE)
					text(x = length(results.to.plot[,p])*1.05, y = sig.effect.levels[pv,p], labels = paste("p =", mark.pvals[pv]), cex = 0.5, adj = 0)
					par(xpd = FALSE)
					}
								
				if(!overlay){
					mtext(pheno.labels[p], cex = 2)
					mtext(paste(pheno.labels[p], "[|Eff|/se]", sep = " "), side = 2, line = 2.5)
					}else{
					mtext("|Eff|/se", side = 2, line = 2.5)	
					}
				axis(2)

				}else{
				
				abline(h = 0)
				if(!overlay){
					mtext(pheno.labels[p], cex = 2)
					mtext(paste(pheno.labels[p], "[Eff]", sep = " "), side = 2, line = 2.5)	
					}else{
					mtext("Eff", side = 2, line = 2.5)		
					}
				
				axis(2)

				# text(1:length(pheno.res), par("usr")[3] - abs((par("usr")[3])*0.1), srt = 90, adj = 1, xpd = TRUE, cex = 0.5)
				} #end case for plotting standardized values
			}
			
			#plot the effect sizes
			for(ch in chr){
				chr.locale <- which(data.obj$chromosome == ch)
				marker.locale <- which(names(pheno.res) %in% chr.locale)
				points(marker.locale, pheno.res[marker.locale], type = plot.type, col = col.mat[,p], pch = 16)
				}
								
		
		if(show.peak.markers){
			points(peak.height, pch = 16, col = "red", cex = 0.4)
			}


		if(show.linkage.blocks){
			null <- lapply(blocks, function(bl) segments(x0 = min(which(names(pheno.res) %in% bl)), y0 = 0, y1 = max(all.vals), lwd = block.delim.lwd, col = block.delim.col, lty = block.delim.lty))
			# par(xpd = TRUE)
			# null <- lapply(blocks, function(bl) points(x = which(names(pheno.res) %in% bl), y = rep((min(all.vals)-(max(all.vals)*0.01)), length(bl)), type = "l", lwd = 3, col = "blue"))
			# single.marker.blocks.locale <- which(lapply(blocks, length) == 1)
			# if(length(single.marker.blocks.locale) > 0){
				# for(i in 1:length(single.marker.blocks.locale)){
					# points(x = which(names(pheno.res) == blocks[[single.marker.blocks.locale[i]]]), y = min(all.vals)-(max(all.vals)*0.01), type = "p", pch = 16, col = "blue", cex = 0.7)
					# }
				# }

			if(mark.block.boundaries){
				null <- lapply(blocks, function(bl) points(x = min(which(names(pheno.res) %in% bl)), y = (min(all.vals)-(max(all.vals)*0.01)), pch = "|", col = "blue", cex = 0.8))
				null <- lapply(blocks, function(bl) points(x = max(which(names(pheno.res) %in% bl)), y = (min(all.vals)-(max(all.vals)*0.01)), pch = "|", col = "blue", cex = 0.8))
				}
			}
			par(xpd = FALSE)

		}
		
				
		if(overlay){
			par(xpd = TRUE)
				legend(x = (length(results.to.plot[,p])*0.94), y = max(all.vals)*1.15, legend = traits, pch = 16, col = trait.colors[1:length(traits)], cex = 0.7)
			par(xpd = FALSE)
		}
		
			

}
