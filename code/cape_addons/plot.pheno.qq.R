

plot.pheno.qq <- function(data.obj, pdf.label = "Phenotype.QQ.plots.pdf"){
	pdf(pdf.label)
	layout(matrix(c(1:4), ncol = 2, byrow = TRUE))
	pheno.pairs <- pair.matrix(colnames(data.obj$pheno))
	for(i in 1:length(pheno.pairs[,1])){
		qqplot(data.obj$pheno[,pheno.pairs[i,1]], data.obj$pheno[,pheno.pairs[i,2]], xlab = pheno.pairs[i,1], ylab = pheno.pairs[i,2])
	}
	dev.off()
	
	}
