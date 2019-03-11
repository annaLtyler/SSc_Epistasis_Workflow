#This function finds the peak marker on the given chromosome(s)
#and returns all markers within the specified effect size drop

peak.markers <- function(data.obj, chr = NULL, effect.drop = 1, pop.out.to.nearest.marker = FALSE){
	
	if(is.null(chr)){
		chr <- unique(data.obj$chromosome[which(data.obj$chromosome != "0")])
		}
	
	all.markers <- 1:length(data.obj$marker.names)
	imp.marker.locale <- grep("loc", data.obj$marker.names)
	if(length(imp.marker.locale) > 0){
		true.marker.locale <- all.markers[-imp.marker.locale]
		}
	
	
	var.to.pheno <- data.obj$max.var.to.pheno.influence
	
	#direction is either 1 or -1 depending on whether you want to go up
	#or down from the starting marker
	find.nearby.markers <- function(chr.markers, max.locale, min.effect, direction){

		all.effect <- chr.markers[max.locale,"|t.stat|"]
		end.reached = 0
		step <- direction
		while(length(which(all.effect <= min.effect)) == 0 && !end.reached){
			test.locale <- max.locale:(max.locale+step)
			if(min(test.locale) < 1 || max(test.locale) > dim(chr.markers)[1]){
				end.reached = 1
				next()
				}				
			all.effect <- chr.markers[test.locale,"|t.stat|"]
			step <- step + direction
			}
		all.markers <- chr.markers[test.locale[which(all.effect > min.effect)],1]
		return(all.markers)		
		}


	all.peaks <- NULL

	for(c in 1:length(chr)){			
		chr.markers <- colnames(data.obj$geno)[which(data.obj$chromosome %in% chr[c])]
		marker.locale <- lapply(var.to.pheno, function(x) which(x[,1] %in% chr.markers))

		for(i in 1:length(var.to.pheno)){
			chr.markers <- var.to.pheno[[i]][marker.locale[[i]],]
			chr.markers <- chr.markers[order(chr.markers[,1]),] #order the markers as they appear on the genome
			max.effect <- max(chr.markers[,"|t.stat|"])
			min.effect <- max.effect[1] - effect.drop
			max.locale <- which(chr.markers[,"|t.stat|"] == max.effect)

			#from this location ratchet out in each direction
			#independently until we are at or below the effect 
			#drop
			upstream <- find.nearby.markers(chr.markers, max.locale, min.effect, direction = 1)
			downstream <- find.nearby.markers(chr.markers, max.locale, min.effect, direction = -1)
			
			all.markers <- sort(unique(c(upstream, downstream)))
			marker.names <- data.obj$marker.names[which(colnames(data.obj$geno) %in% all.markers)]
			
			if(pop.out.to.nearest.marker){
				marker.name.locale <- which(data.obj$marker.names %in% marker.names)
				if(length(which(imp.marker.locale == marker.name.locale[1])) > 0){
					nearest.marker.down <- max(true.marker.locale[which(true.marker.locale < marker.name.locale[1])])
					all.markers <- nearest.marker.down:max(all.markers)
					}
				if(length(which(imp.marker.locale == marker.name.locale[length(marker.name.locale)])) > 0){
					nearest.marker.up <- min(true.marker.locale[which(true.marker.locale > marker.name.locale[length(marker.name.locale)])])
					all.markers <- min(all.markers):nearest.marker.up
					}					
				}
			
			marker.names <- data.obj$marker.names[which(colnames(data.obj$geno) %in% all.markers)]
			marker.chr <- data.obj$chromosome[which(colnames(data.obj$geno) %in% all.markers)]
			marker.pos <- data.obj$marker.location[which(colnames(data.obj$geno) %in% all.markers)]
			phenotype <- rep(names(var.to.pheno)[i], length(all.markers))
			marker.stats <- chr.markers[which(chr.markers[,1] %in% all.markers),((dim(chr.markers)[2]-2):(dim(chr.markers)[2]))]
			
			peak.table <- matrix(c(marker.names, marker.chr, marker.pos, phenotype, marker.stats), ncol = 7, byrow = FALSE)
			all.peaks <- rbind(all.peaks, peak.table)
			}
		}
		
		if(!is.null(all.peaks)){
			colnames(all.peaks) <- c("marker", "marker.chr", "marker.pos", "phenotype", colnames(chr.markers)[((dim(chr.markers)[2]-2):(dim(chr.markers)[2]))])
			write.table(all.peaks, paste("Peak.Markers.", effect.drop, "drop.csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)
			
			pdf(paste("Peak.Markers.", effect.drop, "drop.pdf", sep = ""), width = 12, height = 4*length(var.to.pheno))
			par(mfrow = c(length(var.to.pheno), 1))
			for(i in 1:length(var.to.pheno)){
				pheno.locale <- which(all.peaks[,"phenotype"] == names(var.to.pheno)[i])
				peak.markers <- unique(all.peaks[pheno.locale,"marker"])
				marker.num <- colnames(data.obj$geno)[which(data.obj$marker.names %in% peak.markers)]
				
				effect.table <- var.to.pheno[[i]]
				effect.table <- effect.table[order(effect.table[,1]),]
				marker.locale <- which(effect.table[,1] %in% marker.num)
				plot(effect.table[,"|t.stat|"], type = "h", ylab = names(var.to.pheno)[i])
				points(marker.locale, (effect.table[marker.locale,"|t.stat|"])*1.02, pch = "*", col = "red")
				}
				
			dev.off()		
		}
		
		data.obj$peak.markers <- all.peaks
		return(data.obj)
	
	
}