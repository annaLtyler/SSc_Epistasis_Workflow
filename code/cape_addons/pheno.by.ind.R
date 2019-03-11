#This script plots phenotypes by individual and
#highlights outliers


pheno.by.ind <- function(data.obj, phenotypes = NULL, pdf.label = "Phenotype.by.Ind.pdf"){
	
	
	if(is.null(phenotypes)){
		phenotypes <- colnames(data.obj$pheno)
		}
	
	num.pheno <- length(phenotypes)
	if(num.pheno <= 16){
		layout.mat <- get.layout.mat(num.pheno, "upright")
		}else{
		layout.mat <- get.layout.mat(16, "upright")
		}
	
	plot.pheno <- function(pheno, pheno.name){
		
		pheno.mean <- mean(pheno, na.rm = TRUE)
		pheno.sd <- sd(pheno, na.rm = TRUE)
		
		plot(pheno, main = pheno.name, ylab = pheno.name, xlab = "Individual")
		abline(h = pheno.mean, col = "red", lwd = 3)
		abline(h = (pheno.mean+(2*pheno.sd)), lty = 2, col = "red", lwd = 2)
		abline(h = (pheno.mean-(2*pheno.sd)), lty = 2, col = "red", lwd = 2)

		}

	pdf(pdf.label, height = dim(layout.mat)[1]*3, width = dim(layout.mat)[2]*3)
	layout(layout.mat)
	for(p in phenotypes){
		plot.pheno(data.obj$pheno[,p], p)
		}
	dev.off()
	
	
	
}