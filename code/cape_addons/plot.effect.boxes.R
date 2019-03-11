#This function plots significant interactions
#bin.geno can be set to something like c(0, 0.5, 1)
#if you have continuous genotypes and want to bin
#them into something recognizable.

plot.effect.boxes <- function(data.obj, p.or.q = 0.05, pdf.label = "Interaction.BoxPlots.pdf", raw.pheno.normalized = TRUE, bin.geno = NULL, by.blocks = FALSE, r2.thresh = 0.8){
	
	var.influences <- data.obj$var.to.var.p.val
	
	pheno.names <- colnames(data.obj$pheno)
	num.pheno <- length(pheno.names)
	num.et <- dim(data.obj$ET)[2]
	et.names <- colnames(data.obj$ET)
	
	geno <- data.obj$geno

	if(by.blocks){
		old.r2.thresh <- data.obj$r2.thresh
		if(is.null(old.r2.thresh) || old.r2.thresh != r2.thresh){
			data.obj <- linkage.blocks(data.obj, p.or.q = p.or.q, r2.thresh = r2.thresh)
			}
		blocks <- data.obj$linkage.blocks.collapsed
				
		new.geno <- sapply(blocks, function(x) rowMeans(geno[,x,drop=FALSE]))
		geno <- new.geno
		colnames(geno) <- sig.markers <- 1:length(blocks)

		if(!is.null(bin.geno)){
			geno <- t(apply(geno, 1, function(x) bin.vector(x, bins = bin.geno)))
			colnames(geno) <- 1:length(blocks)
			}
		
		marker.names <- names(blocks)
		all.pairs <- pair.matrix(1:length(blocks))
		sig.int.locale <- 1:length(blocks)
		
		}else{
		
		if(!is.null(bin.geno)){
			geno <- t(apply(geno, 1, function(x) bin.vector(x, bins = bin.geno)))
			colnames(geno) <- colnames(data.obj$geno)
			}
		marker.names <- data.obj$marker.names
		all.pairs <- var.influences[,1:2]
		var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))
		sig.int.locale <- which(var.influences[,var.sig.col] <= p.or.q)
		}
			


	if(raw.pheno.normalized){
		all.pheno <- cbind(data.obj$pheno, data.obj$ET)
		}else{
		#use the uncentered, unnormalized data
		all.pheno <- cbind(data.obj$raw.pheno, data.obj$ET)
		}
		
	all.pheno.names <- c(pheno.names, et.names)
	total.response <- num.et + num.pheno
	
	if(is.null(var.influences)){
		stop("variant-to-variant influences must be calculated before plotting.")
		}

	
	if(length(sig.int.locale) == 0){
		
		pdf(pdf.label)
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(0,1))
		text(0.5, 0.5, "No Significant Interactions")
		dev.off()
		}else{
		
	
		pdf(pdf.label, width = 12, height = total.response*3)
		layout(matrix(1:(total.response*3), nrow = total.response,  ncol = 3, byrow = TRUE), widths = c(1,1,2))
		for(i in sig.int.locale){
			marker1 <- geno[,as.character(all.pairs[i,1])]
			marker2 <- geno[,as.character(all.pairs[i,2])]
			marker.table <- cbind(marker1, marker2)
			marker.labels <- apply(marker.table, 1, function(x) paste(x, collapse = ","))
			not.na <- intersect(which(!is.na(marker1)), which(!is.na(marker2)))
			for(ph in 1:total.response){
				boxplot(all.pheno[,ph]~marker1, notch = TRUE, ylab = pheno.names[ph], main = paste(marker.names[all.pairs[i,1]], ", ", all.pheno.names[ph], sep = ""))	
				boxplot(all.pheno[,ph]~marker2, notch = TRUE, ylab = pheno.names[ph], main = paste(marker.names[all.pairs[i,2]], ", ", all.pheno.names[ph], sep = ""))	
				boxplot(all.pheno[not.na,ph]~marker.labels[not.na], notch = TRUE, ylab = pheno.names[ph], main = paste(marker.names[all.pairs[i,1]], ", ", marker.names[all.pairs[i,2]], "\n", all.pheno.names[ph], sep = ""))	
				}
						
			}

		dev.off()
		}
}

