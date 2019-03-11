

prop.alt.vs.pheno.var <- function(data.obj, markers = NULL, phenotype = NULL, covar = NULL, plot.title = NULL){
	
	
	if(is.null(markers)){
		markers <- data.obj$marker.names
		}
	
	pheno <- data.obj$pheno
	
	if(!is.null(phenotype)){
		pheno.locale <- match(phenotype, colnames(pheno))
		pheno <- pheno[,pheno.locale,drop=FALSE]
		}

	#pull out the effects of the covariate(s)
	if(!is.null(covar)){
		covar.locale <- match(covar, data.obj$marker.names)
		covar.mat <- data.obj$geno[,covar.locale, drop = FALSE]
		#mean center the covariates
		centered <- apply(covar.mat, 2, function(x) x - mean(x, na.rm = TRUE))
		models <- apply(pheno, 2, function(x) lm(x~covar.mat))
		resids <- lapply(models, residuals)
	
		na.locale <- unique(which(is.na(cbind(pheno, covar.mat)), arr.ind = TRUE)[,1])
		not.na <- c(1:dim(pheno)[1])
		if(length(na.locale) > 0){
			not.na <- not.na[-na.locale]
			}
		for(ph in 1:dim(pheno)[2]){
			pheno[not.na,ph] <- resids[[ph]]
			}
		}
	
	marker.locale <- match(markers, data.obj$marker.names)
	just.wanted.markers <- data.obj$geno[,marker.locale]
	
	#find the amount of the alternate allele each individual has at these markers
	total.alternate <- apply(just.wanted.markers, 1, function(x) sum(x)/length(x))
	
	#plot the phenotype vs percent alternate allele
	layout.mat <- get.layout.mat(dim(pheno)[2]*2, "landscape")
	
	# quartz(width = dim(layout.mat)[2]*4, height = dim(layout.mat)[1]*4)
	
	layout(layout.mat)
	for(i in 1:dim(pheno)[2]){
		model <- lm(pheno[,i]~total.alternate)
		r.squared <- summary(model)$r.squared
		pval <- summary(model)$coefficients[2,"Pr(>|t|)"]

		plot(total.alternate, pheno[,i], main = paste("Phenotype vs. Percent C3H\n", colnames(pheno)[i], "\nR2 =", signif(r.squared, 2), ", p =", signif(pval, 2), sep = ""), xlab = "Percent C3H", ylab = colnames(pheno)[i])
		abline(model, col = "red", lwd = 3)
		}
	

	#plot histograms along a sliding window
	bins <- seq(0,1,0.05)
	
	# quartz(width = dim(layout.mat)[2]*4, height = dim(layout.mat)[1]*4)
	# layout(layout.mat)
	for(ph in 1:dim(pheno)[2]){
		plot.new()
		plot.window(xlim = c(1,length(bins)), ylim = c(-4,4))
		for(i in 1:(length(bins)-2)){
			ind.locale <- intersect(which(total.alternate > bins[i]), which(total.alternate <= bins[(i+2)]))
			if(length(ind.locale) > 1){
				mean.pheno <- mean(pheno[ind.locale,ph], na.rm = TRUE)
				to.plot <- pheno[ind.locale,ph] - mean.pheno
				boxplot(to.plot, at = i, add = TRUE, xlab = colnames(pheno)[i])
				}
			}
		mtext("Binned, Mean-Centered Phenotype\n vs. Percent C3H", side = 3, cex = 0.7)
		mtext(colnames(pheno)[ph], side = 1, line = 1.5)
		}
	
	
	mtext(plot.title, outer = TRUE, line = -1)

	# #also look at percent heterozygosity
	# binned.geno <- apply(just.wanted.markers, 2, function(x) bin.vector(x, bins = bin.geno))
	# percent.het <- apply(binned.geno, 1, function(x) length(which(x == 0.5))/length(x))
	
	
	# for(ph in 1:dim(pheno)[2]){
		# plot(percent.het, pheno[,ph])
		# }

	}