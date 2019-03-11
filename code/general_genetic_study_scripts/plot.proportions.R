plot.proportions <- function(data.obj, geno.obj, snp1, snp2, trait, allele.which = 2, cex.axis = 1, cex.names = 1, plot.type = c("bar", "line"), trait.type = c("ET", "Normalized", "Raw"), covar = NULL){
	
	geno <- get.geno(data.obj, geno.obj)
	cor.pheno <- get.pheno(data.obj, trait.type[1], covar = covar)
		
	snp1.locale <- which(dimnames(geno)[[3]] == snp1)
	snp2.locale <- which(dimnames(geno)[[3]] == snp2)
	trait.locale <- which(colnames(data.obj$pheno) == trait)
	if(length(trait.locale) == 0){stop("Can't find the trait:", trait)}
	pheno <- cor.pheno[,trait.locale]
	
	has.snp1 <- which(geno[,allele.which,snp1.locale] > 0)
	no.snp1 <- which(geno[,allele.which,snp1.locale] == 0)
	
	has.snp2 <- which(geno[,allele.which,snp2.locale] > 0)
	no.snp2 <- which(geno[,allele.which,snp2.locale] == 0)
	
	no.exposure <- intersect(no.snp1, no.snp2)
	both.exposure <- intersect(has.snp1, has.snp2)
	means <- c(mean(pheno[no.exposure]), mean(pheno[has.snp1]), mean(pheno[has.snp2]), mean(pheno[both.exposure]))
	n <- c(length(no.exposure), length(has.snp1), length(has.snp2), length(both.exposure))
	
	means <- means - means[1]
	plot.height <- max(means) - min(means)
	max.mean <- max(means)
	bump.up <- plot.height*0.05
	bump.down <- plot.height*0.13
	bump.left <- 0
	par(mar = c(8,8,4,2))
	
	plot.min <- min(means)*1.2
	plot.max <- max(means)*1.2
	
	
	if(plot.type == "bar"){
		a <- barplot(means, main = trait, ylim = c(plot.min, plot.max), cex.axis = cex.axis)
		}else{
		plot(means, type = "b", ylim = c(plot.min, plot.max), axes = FALSE, ylab = "", xlab = "")
		axis(2)
		a <- matrix(c(1:4), ncol = 1)
		}
	par(xpd = TRUE)
	#write the SNP names and their presence/absence markers
	text(x = c(bump.left, a[,1]), y = plot.min-bump.down, labels = c(snp1, 0, 1, 0, 1), cex = cex.names, adj  = c(1,0.5,0.5,0.5))
	text(x = c(bump.left, a[,1]), y = plot.min-(bump.down*2), labels = c(snp2, 0, 0, 1, 1), cex = cex.names, adj  = c(1,0.5,0.5,0.5))
	# text(x = a[,1], y = means+bump.up, labels = signif(means, 2), cex = cex.names)
	text(x = a[,1], y = rep(plot.max+(bump.up*2), length(means)), labels = n, cex = cex.names)
	par(xpd = FALSE)
	
	
	
	
}