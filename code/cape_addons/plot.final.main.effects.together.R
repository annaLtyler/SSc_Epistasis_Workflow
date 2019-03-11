#This script allows plotting main effects from different crosses
#on the same set of axes as long as the crosses are comparable
#the multiple cross object is generated

plot.final.main.effects.together <- function(multiple.cross.obj, chr = NULL, traits = NULL, standardized = TRUE, mark.chr = TRUE, plot.type = "h", cross.colors = NULL, mark.pvals = c(0.01, 0.05)){
	
	num.crosses <- length(multiple.cross.obj) - 3	
	#order the results by marker order
	D1.results <- vector("list", num.crosses)
	for(i in 1:num.crosses){
		D1.results[[i]] <- multiple.cross.obj[[i]]
		D1.results[[i]] <- lapply(D1.results[[i]], function(x) x[order(x[,1]),])
		}
	
	marker.names <- multiple.cross.obj$marker.names
	all.marker.num <- 1:length(marker.names)

	if(is.null(chr)){
		all.chr <- NULL
		for(i in 1:length(D1.results)){
			all.chr <- c(all.chr, unique(multiple.cross.obj$chromosome[which(all.marker.num %in% D1.results[[i]][[1]][,1])]))
			}
		all.chr <- sort(unique(all.chr))
		}

	if(is.null(traits)){
		all.traits <- NULL
		for(i in 1:length(D1.results)){
			all.traits <- c(all.traits, names(D1.results[[i]]))
			}
		all.traits <- unique(all.traits)	
		}
		
	#get the colors of the points to use for each cross
	col.list <- vector("list", length = num.crosses)

		if(is.null(cross.colors)){
			cross.colors <- c("black", "blue", "purple", "darkgreen")
			}
		if(length(cross.colors) < length(D1.results)){
		 	cross.colors <- rep(cross.colors, length(D1.results)/4)
		 	cross.colors <- cross.colors[1:length(D1.results)]
			}
		for(i in 1:length(D1.results)){
			col.list[[i]] <- rep(cross.colors[i], dim(D1.results[[i]][[1]])[1])
			}


	chr.locale <- which(multiple.cross.obj$chromosome %in% all.chr)
	marker.chr <- lapply(D1.results, function(x) multiple.cross.obj$chromosome[which(all.marker.num %in% x[[1]][,1])])
	markers.used <- lapply(D1.results, function(x) which(all.marker.num %in% x[[1]][,1]))
	markers.which <- lapply(markers.used, function(x) intersect(chr.locale, x))
	results.rows <- vector("list", length(markers.which))
	for(i in 1:length(markers.which)){
		results.rows[[i]] <- which(D1.results[[i]][[1]][,1] %in% all.marker.num[markers.which[[i]]])
		}
	results.el <- lapply(D1.results, function(x) which(names(x) %in% all.traits))
	results.to.plot <- vector("list", length(D1.results)); names(results.to.plot) <- cross.names
	for(i in 1:length(D1.results)){
		for(r in results.el[[i]]){
			if(standardized){
				results.to.plot[[i]] <- cbind(results.to.plot[[i]], D1.results[[i]][[r]][results.rows[[i]],"|t.stat|"])
				}else{
				results.to.plot[[i]] <- cbind(results.to.plot[[i]], D1.results[[i]][[r]][results.rows[[i]],"coef"])
				}	
			}
		}

	for(i in 1:length(results.to.plot)){
		rownames(results.to.plot[[i]]) <- D1.results[[i]][[1]][results.rows[[i]],1]
		colnames(results.to.plot[[i]]) <- names(D1.results[[i]])[results.el[[i]]]
		}

	#find the effect sizes where the specified corrected p values lie
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(D1.results[[1]][[1]]))))

		sig.effect.levels <- vector("list", length(D1.results))
		for(i in 1:length(D1.results)){
			for(pv in 1:length(mark.pvals)){
				if(standardized){
					sig.effect.levels[[i]] <- rbind(sig.effect.levels[[i]], unlist(lapply(D1.results[[i]], function(x) min(x[which(x[,var.sig.col] <= mark.pvals[pv]),"|t.stat|"]))))
					}else{
					sig.effect.levels[[i]] <- rbind(sig.effect.levels[[i]], unlist(lapply(D1.results[[i]], function(x) min(x[which(x[,var.sig.col] <= mark.pvals[pv]),"coef"]))))	
					}
				}
			}


	layout.mat <- matrix(1:length(all.traits), ncol = 1)

	layout(layout.mat)	
	
	for(p in 1:length(all.traits)){
		# dev.new(width = plot.width, height = plot.height)
		pheno.res <- lapply(results.to.plot, function(x) x[,p])
		
		if(standardized){
			all.vals <- c(unlist(results.to.plot), 0)	
			}else{
			all.vals <- unlist(results.to.plot)	
			}
			

		#create a new window for each trait
		par(mar = c(3, 4, 3, 2) + 0.1)
		plot.new()
		plot.window(xlim = c(0, length(multiple.cross.obj$chromosome)), ylim = c(min(all.vals), max(all.vals[is.finite(all.vals)])))
		
			#shade the chromosome regions based on all the chromosomes
			if(mark.chr){
				par(xpd = TRUE)
				for(ch in 1:length(all.chr)){
					x.min <- min(which(multiple.cross.obj$chromosome == all.chr[ch])); x.max <- max(which(multiple.cross.obj$chromosome == all.chr[ch]))
					if(ch %% 2 == 0){
						polygon(x = c(x.min, x.min, x.max, x.max), y = c(min(all.vals), max(all.vals), max(all.vals), min(all.vals)), col = "lightgray", border = NA)
						}
					if(all.chr[ch] == 0){
						text(x = x.max, y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = "Cov.", cex = 0.5, adj = 0)
						}else{
						text(x = mean(c(x.min, x.max)), y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = all.chr[ch], cex = 0.5)
						}
					}
				par(xpd = FALSE)
				}
		
				abline(h = 0)




			if(standardized){

				#add lines to indicate the significance thresholds
				for(pv in 1:length(mark.pvals)){
					for(i in 1:length(sig.effect.levels)){
						abline(h = sig.effect.levels[[i]][pv,p], lty = pv, col = col.mat[1,i])
						par(xpd = TRUE)
						# text(x = length(results.to.plot[[i]][,p])*1.05, y = sig.effect.levels[[i]][pv,p], labels = paste(cross.names[i], "p =", mark.pvals[pv]), cex = 0.5, adj = 0)
						par(xpd = FALSE)
						}
					}
								
					mtext(all.traits[p], cex = 2)
					mtext(paste(colnames(results.to.plot)[p], "[|Eff|/se]", sep = " "), side = 2, line = 2.5)
					axis(2)

				}else{
				
				abline(h = 0)

				mtext(all.traits[p], cex = 2)
				mtext(paste(all.traits[p], "[Eff]", sep = " "), side = 2, line = 2.5)	
				
				axis(2)

				
			}
			
			#plot the effect sizes from each cross
			for(i in 1:length(pheno.res)){
				for(ch in all.chr){
					chr.locale <- which(marker.chr[[i]] == ch)
					if(length(chr.locale) > 0){
						points(x = markers.which[[i]][chr.locale], pheno.res[[i]][chr.locale], type = plot.type, col = col.list[[i]], pch = 16)
						}
					}
				}
		plot.dim <- par("usr")
		par(xpd = TRUE)
		legend(x = plot.dim[2]-(plot.dim[2]-plot.dim[1])*0.15, y = plot.dim[4]*1.15, legend = cross.names, pch = 16, col = cross.colors[1:length(cross.names)], cex = 0.7)
		par(xpd = FALSE)
					

		}
		
				


		
			

}
