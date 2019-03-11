#This function plots the raw value of a continuous covariate
#by the genotype of an interacting marker.
#continuous.covar <- "Final.IGF.1"; interacting.marker = c("Chr7_loc10", "Chr7_loc20", "Chr7_D7Mit228", "Chr7_loc25", "Chr9_loc20", "Chr9_D9Mit259", "Chr9_loc25", "Chr9_loc30", "Chr9_D9Mit12", "Chr9_loc35", "Chr9_loc40")
#plot.continuous.interaction(cross, continuous.covar = continuous.covar, interacting.marker = interacting.marker)

plot.continuous.interaction <- function(data.obj, continuous.covar, interacting.marker = NULL, rebin.geno = TRUE, genotype.bins = c(0, 0.5, 1), genotype.col = c("blue", "purple", "forestgreen"), show.linear.fit = TRUE, show.spline.fit = FALSE){
	
	
	covar.locale <- which(data.obj$marker.names %in% continuous.covar)
	if(length(covar.locale) != length(continuous.covar)){
		found <- data.obj$marker.names[covar.locale]
		not.found <- setdiff(continuous.covar, found)
		cat("Please check spelling. I couldn't find the following covariates:\n")
		cat(not.found, sep = "\n")
		return(NULL)
		}
	
	
	if(is.null(interacting.marker)){
		interacting.marker <- data.obj$marker.names[-covar.locale]
		}
		
	marker.locale <- which(data.obj$marker.names %in% interacting.marker)
	if(length(marker.locale) != length(interacting.marker)){
		found <- data.obj$marker.names[marker.locale]
		not.found <- setdiff(interacting.marker, found)
		cat("Please check spelling. I couldn't find the following markers:\n")
		cat(not.found, sep = "\n")
		return(NULL)
		}
	
	num.pheno <- dim(data.obj$pheno)[2]
	
	layout.matrix <- matrix(1:((num.pheno)*2 + (num.pheno*length(genotype.bins))), nrow = num.pheno, byrow = TRUE)
	
	pdf(paste("Interactions.No.Binning.", interacting.marker, "_", continuous.covar, ".pdf", sep = ""), width = dim(layout.matrix)[2]*5, height = dim(layout.matrix)[1]*4)
	layout(layout.matrix)
	# layout.show(16)
	for(i in 1:length(covar.locale)){
		for(j in 1:length(marker.locale)){

			all.geno <- data.obj$geno[,marker.locale[j]]
			if(rebin.geno){
				binned.geno <- bin.vector(all.geno, bins = genotype.bins)
				}else{
				binned.geno <- all.geno
				}
			pt.col <- rep(NA, length(binned.geno))
			for(g in 1:length(genotype.bins)){
				pt.col[which(binned.geno == genotype.bins[g])] <- genotype.col[g]
				}
			for(p in 1:dim(data.obj$pheno)[2]){
				par(mar = c(4,4,3,2))
				plot(data.obj$geno[,covar.locale[i]], data.obj$pheno[,p], col = pt.col, xlab = continuous.covar[i], ylab = colnames(data.obj$pheno)[p], main = paste(colnames(data.obj$pheno)[p], data.obj$marker.names[marker.locale[j]], sep = "\n"))
				
				model.list <- vector(mode = "list", length = length(genotype.bins))
				
				for(cl in 1:length(genotype.bins)){
					cl.locale <- which(pt.col == genotype.col[cl])
					plot(data.obj$geno[cl.locale,covar.locale[i]], data.obj$pheno[cl.locale,p], col = genotype.col[cl], xlab = continuous.covar[i], ylab = colnames(data.obj$pheno)[p], main = paste("genotype:", genotype.bins[cl]), ylim = c(min(data.obj$pheno[,p], na.rm = TRUE), max(data.obj$pheno[,p], na.rm = TRUE)))
					if(show.spline.fit){
						val.mat <- cbind(data.obj$geno[cl.locale,covar.locale[i]], data.obj$pheno[cl.locale,p])
						good.vals <- val.mat[which(!is.na(rowSums(val.mat))),]
						test <- smooth.spline(good.vals)
						points(test, type = "l", lty = 2, lwd = 2)
						}
					if(show.linear.fit){
						model <- lm(data.obj$pheno[cl.locale,p]~data.obj$geno[cl.locale,covar.locale[i]])
						pval <- anova(model)$"Pr(>F)"[1]
						legend("topright", legend = paste("p =", signif(pval, 3)), lty = 1, lwd = 2)
						abline(model)
						model.list[[cl]] <- model
						}
					}

				plot.new()
				plot.window(xlim = c(0,1), ylim = c(-1,1))
				axis(1, cex.axis = 2); axis(2, cex.axis = 2)
				for(cl in 1:length(genotype.bins)){
					abline(model.list[[cl]], lty = c(length(genotype.bins):1)[cl], lwd = 3)
					}
				legend("topleft", lty = c(length(genotype.bins):1), legend = genotype.bins, cex = 2, lwd = 3)
				mtext(text = "Slope Comparison")



				# par(xpd = TRUE)
				# legend("topright", legend = genotype.bins, inset = c(-0.18, 0), col = genotype.col, pch = 16)
				# par(xpd = FALSE)
				}
			}		
		
		}
	
	dev.off()

	
}