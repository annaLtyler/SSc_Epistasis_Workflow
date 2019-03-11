
plot.pheno.dist <- function(data.obj, pdf.label = "Phenotype.Distributions.pdf"){

	pdf(pdf.label)
	layout(matrix(c(1:4), ncol = 2, byrow = TRUE))
	for(i in 1:dim(data.obj$pheno)[2]){
		hist(data.obj$pheno[,i], main = colnames(data.obj$pheno)[i], xlab = colnames(data.obj$pheno)[i])
		}
	dev.off()

	}