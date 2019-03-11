#This function plots all context-dependent effects of a marker
#the argument color.over.perc.slope determines the percentage 
#of the maximum slope that constitutes an actual change. Lines 
#with less than this slope will be colored gray

plot.main.effects.all.context <- function(data.obj, p.or.q = 0.05, color.over.perc.slope = 1, pdf.label = "Main.Effects.All.Pair.Contexts.pdf", standardized = TRUE, ylim = NULL, bin.geno = c(0,0.5,1), by.blocks = TRUE, r2.thresh = 0.5, alt.bin.var = NULL, alt.bins = NULL){

	
	var.influences <- data.obj$max.var.to.pheno.influence
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences[[1]]))))
	sig.data <- lapply(var.influences, function(x) x[which(as.numeric(x[,var.sig.col]) <= p.or.q),,drop = FALSE])
	sig.markers <- sort(unique(unlist(lapply(sig.data, function(x) x[,1]))))

	phenotypes <- names(var.influences)	
	num.pheno <- length(phenotypes)
	geno <- data.obj$geno
		
	if(is.null(bin.geno)){
		xlim <- c(0, max(geno))
		}else{
		xlim <- c(0, max(bin.geno))	
		}

	if(by.blocks){
		old.r2.thresh <- data.obj$r2.thresh
		if(is.null(old.r2.thresh) || old.r2.thresh != r2.thresh){
			cat("Recalculating linkage blocks...\n")
			data.obj <- linkage.blocks(data.obj, p.or.q = p.or.q, r2.thresh = r2.thresh)
			}
		blocks <- data.obj$linkage.blocks.collapsed
				
		new.geno <- sapply(blocks, function(x) rowMeans(geno[,x,drop=FALSE]))
		geno <- new.geno
		colnames(geno) <- sig.markers <- 1:length(blocks)

		#bin the genotypes if this has been specified
		if(!is.null(bin.geno)){
			if(!is.null(alt.bin.var)){
				#but leave out the columns for alternate binning
				alt.bin.locale <- match(alt.bin.var, names(blocks))
				alt.bin.vals <- matrix(geno[,alt.bin.locale], ncol = length(alt.bin.var), byrow = FALSE)
				}
			geno <- t(apply(geno, 1, function(x) bin.vector(x, bins = bin.geno)))
			if(!is.null(alt.bin.var)){
				geno[,alt.bin.locale] <- apply(alt.bin.vals, 2, function(x) bin.vector(x, bins = alt.bins))
				}
			colnames(geno) <- 1:length(blocks)
			}
		
		marker.names <- names(blocks)
		all.pairs <- pair.matrix(1:length(blocks))

		}else{
		
		if(!is.null(bin.geno)){
			if(!is.null(alt.bin.var)){
				alt.bin.locale <- match(alt.bin.var, data.obj$marker.names)
				alt.bin.vals <- geno[,alt.bin.locale]
				}
			geno <- t(apply(geno, 1, function(x) bin.vector(x, bins = bin.geno)))
			if(!is.null(alt.bin.var)){
				geno[,alt.bin.locale] <- alt.bin.vals
				}
			colnames(geno) <- colnames(data.obj$geno)
			}
		marker.names <- data.obj$marker.names
		all.pairs <- data.obj$var.to.var.p.val[,1:2]
	
		}



	plot.effects <- function(marker.num){

		marker.name <- marker.names[marker.num]
		
		#find all other markers this one gets paired with
		paired.locale <- c(which(all.pairs[,1] == marker.num), which(all.pairs[,2] == marker.num))
		paired.markers <- unique(as.vector(all.pairs[paired.locale,]))
		paired.markers <- sort(paired.markers[which(paired.markers != marker.num)])
	
		if(length(paired.markers) == 0){
			return("No markers were paired with this one")
			}			
		
		#go through each of these pairings, condition on the second marker
		#and show the effect of the marker of interest on the phenotype

		marker.geno <- unique(geno[,as.character(marker.num)]) #get all unique genotypes of the marker we are looking at
		marker.geno <- sort(marker.geno[which(!is.na(marker.geno))])
		
		#find all locations where the marker we are looking at takes each of its genotype values
		geno.locale <- apply(matrix(marker.geno, ncol = 1), 1, function(x) which(geno[,as.character(marker.num)] == x))

		#if the genotypes show up in equal numbers, the value returned will be a matrix.
		#make sure it's a list
		if(class(geno.locale) == "matrix"){
			geno.list <- vector(mode = "list", length = dim(geno.locale)[2])
			for(g in 1:length(marker.geno)){
				geno.list[[g]] <- geno.locale[,g]
				}
			geno.locale <- geno.list
			}
		names(geno.locale) <- marker.geno
		
		#for each phenotype make a table of the main effect of
		#the marker of interest in each genetic context		
		for(p in 1:num.pheno){
			# cat("\t", phenotypes[ph], "\n")
			effect.mat <- NULL
			
			#collect the data for each marker the marker of interest is paired with
			for(i in 1:length(paired.markers)){
				# report.progress(i, length(paired.markers))
				paired.marker.geno <- geno[,as.character(paired.markers[i])]
				paired.marker.geno <- sort(unique(paired.marker.geno[which(!is.na(paired.marker.geno))]))
				
				#find out where the paired marker takes on each of its unique genotypes
				paired.geno.locale <- apply(matrix(paired.marker.geno, ncol = 1), 1, function(x) which(geno[,as.character(paired.markers[i])] == x))
				names(paired.geno.locale) <- paired.marker.geno
				
				#for each genotype of the paired marker
				#find the effect of the marker of interest
				for(j in 1:length(paired.geno.locale)){
					overlap <- lapply(geno.locale, function(x) intersect(x, paired.geno.locale[[j]]))
					pheno.vals <- lapply(overlap, function(x) data.obj$pheno[x,p])
					pheno.means <- sapply(pheno.vals, function(x) mean(x, na.rm = TRUE))
					pheno.se <- sapply(pheno.vals, function(x) sd(x)/sqrt(length(x)))
					pheno.sd <- sapply(pheno.vals, sd)
					if(standardized){
						effect.mat <- rbind(effect.mat, pheno.means/pheno.sd)
					}else{
						effect.mat <- rbind(effect.mat, pheno.means)
						}
					} #end looping over the paired marker genotypes

				}#end looping over the paired markers
			
			bad.vals <- unique(as.vector(unlist(apply(effect.mat, 2, function(x) which(!is.finite(x))))))
			if(length(bad.vals) > 0){
				effect.mat <- effect.mat[-bad.vals,,drop=FALSE]
				}
			
			#take out duplicate rows
			effect.mat <- unique(effect.mat, MARGIN = 1)
			effect.cols <- rep(NA, dim(effect.mat)[1])
			effect.change <- apply(effect.mat, 1, function(x) x[length(x)] - x[1])
			effect.cols[which(effect.change < 0)] <- "cornflowerblue"
			effect.cols[which(effect.change > 0)] <- "deeppink"
			no.col.val <- max(abs(effect.change))*(color.over.perc.slope/100)
			effect.cols[which(abs(effect.change) < no.col.val)] <- "mediumvioletred"
			
			plot.new()
			if(is.null(ylim)){
				plot.window(xlim = xlim, ylim = c(min(effect.mat, na.rm = TRUE), max(effect.mat, na.rm = TRUE)))
				}else{
				plot.window(xlim = xlim, ylim = ylim)	
				}
			#plot.window(xlim = c(0,1), ylim = c(min(data.obj$pheno[,p], na.rm = TRUE), max(data.obj$pheno[,p], na.rm = TRUE)))
			mtext(paste(phenotypes[p], marker.name, sep = "\n"))
			mtext(phenotypes[p], side = 2, line = 2.5)
			mtext(marker.name, side = 1, line = 2.5)
			axis(2); axis(1, at = as.numeric(marker.geno))
			for(e in 1:length(effect.mat[,1])){
				points(as.numeric(marker.geno), y = effect.mat[e,], type = "l", col = effect.cols[e])				
				}
			if(is.null(ylim)){
				suppressWarnings(boxplot(effect.mat, main = paste(phenotypes[p], marker.name, sep = "\n"), notch = TRUE, xlab = "Genotype"))
				}else{
				suppressWarnings(boxplot(effect.mat, main = paste(phenotypes[p], marker.name, sep = "\n"), notch = TRUE, xlab = "Genotype", ylim = ylim))	
				}
			}
		}
	
	num.markers.per.page = num.pheno
	if(length(sig.markers) <= num.markers.per.page){
		layout.mat <- matrix(1:(length(sig.markers)*(num.pheno*2)), ncol = num.pheno*2, byrow = TRUE)
		}else{
		layout.mat <- matrix(1:(num.markers.per.page*(num.pheno*2)), ncol = num.pheno*2, byrow = TRUE)
		}
	
	pdf(file = pdf.label, width = 4*dim(layout.mat)[2], 4*dim(layout.mat)[1])
	layout(layout.mat)
	for(s in 1:length(sig.markers)){
		#print(s)
		# cat("\nMarker", names(blocks)[s], "\n")
		report.progress(s, length(sig.markers))
		plot.effects(marker.num = sig.markers[s])
		}
	cat("\n")
	dev.off()
	
	
}