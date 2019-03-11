#This is a wrapper for running an R/qtl analog of our cape analysis
# filename = "NON_NZO_Reifsnyder_pgm_CAPE_lett.csv"; singlescan.perm = 500; pairscan.perm = 1000; pheno.col = 1:4; covar = "mom"; data.dir = "~/Documents/Data/Reifsnyder"; results.dir = "./Rqtl_results"

qtl.cape.analog <- function(filename = "NON_NZO_Reifsnyder_pgm_CAPE_lett.csv", singlescan.perm = 500, pairscan.perm = 0, pheno.col = 1:4, covar = "mom", data.dir = ".", results.dir = "./Rqtl_results"){

library(qtl)
setwd(data.dir)

#====================================================================================
#additional functions
#====================================================================================

get.sig.lod <- function(pairscan.results, pairscan.results.perm, model.which, alpha = 0.05){
	phenos <- dimnames(pairscan.results[[1]])[[3]]
	all.sig <- vector(mode = "list", length = length(phenos))
	names(all.sig) <- phenos
	for(ph in 1:length(phenos)){
		ph.results <- summary(pairscan.results, lodcol = ph)
		model.locale <- grep(model.which, colnames(ph.results))
		sig.level <- summary(pairscan.results.perm, alpha = alpha)[[model.which]][ph]
		sig.results <- which(ph.results[,model.locale] >= sig.level)
		if(length(sig.results) > 0){
			all.sig[[ph]] <- ph.results[sig.results,]
			}else{
				all.sig[[ph]] <- paste("No significant pairs for the ", model.which, " model at alpha = ", alpha, ".", sep = "")
			}
		}
	return(all.sig)
}
#====================================================================================


cross <- read.cross(format = "csv", file = filename)

if(!file.exists(results.dir)){
	system(paste("mkdir", results.dir))
	}

setwd(results.dir)

cross <- sim.geno(cross)
singlescan.results <- scanone(cross, pheno.col = pheno.col, method = "imp", addcovar = cross$pheno$covar)
save(singlescan.results, file = "ondD.results.RData")
singlescan.results.perm <- scanone(cross, pheno.col = pheno.col, method = "imp", addcovar = cross$pheno$covar, n.perm = singlescan.perm)
save(singlescan.results.perm, file = "singlescan.results.perm.RData")
perm.summ <- summary(singlescan.results.perm)


pdf("Rqtl.singlescan.pdf", width = 12, height = 12)
	layout(matrix(1:4, nrow = 4))
	for(i in pheno.col){
		plot(singlescan.results, lodcolumn = i, main = colnames(cross$pheno)[i])
		abline(h = perm.summ[1,i], lty = 2)
		abline(h = perm.summ[2,i], lty = 1)
		}
	dev.off()
	
	
pairscan.results <- scantwo(cross, pheno.col = pheno.col, method = "imp", addcovar = cross$pheno$mom)
saveRDS(pairscan.results, file = "pairscan.results.RData")

pdf("pairscan.results.pdf")
for(i in pheno.col){
	plot(pairscan.results, lodcol = i)
	}
dev.off()


if(pairscan.perm > 0){
	pairscan.results.perm <- scantwo(cross, pheno.col = pheno.col, method = "imp", addcovar = cross$pheno$mom, n.perm = pairscan.perm)
	saveRDS(pairscan.results.perm, file = "pairscan.results.perm.RData")
	# summary(pairscan.results.perm, alpha = 0.05)
	pdf("Null.Distributions.pdf")
	plot(pairscan.results.perm)
	dev.off()
	}


std.thresh <- c(9.1, 7.1, 6.3, 6.3, 3.3)

# setwd(results.dir)
# pairscan.results <- readRDS("pairscan.results.RData")
# pairscan.results.perm <- readRDS("pairscan.results.perm.RData")

sig.int <- get.sig.lod(pairscan.results, pairscan.results.perm, model.which = "int", alpha = 0.05)
unlink("QTL.Scantwo.Results.txt")
sink(file = "QTL.Scantwo.Results.txt", append = TRUE)
print(sig.int)
sink()

unlink("QTL.Scantwo.Results.Broman.txt")
sink(file = "QTL.Scantwo.Results.Broman.txt", append = TRUE)
for(i in pheno.col){
	cat(paste("\n\n", colnames(cross$pheno)[i], "\n"))
	# print(summary(pairscan.results, perm = pairscan.results.perm, alphas = rep(0.05, 5), lodcol = i))
	# print(summary(pairscan.results, perm = pairscan.results.perm, alphas = c(0.05, 0.05, 0, 0.05, 0.05), lodcol = i, pvalues = TRUE))
	# print(summary(pairscan.results, perm = pairscan.results.perm, alphas = c(0.05, 0.05, 0.05, 0.05, 0.05), lodcol = i, what = "int"))
	if(pairscan.perm > 0){
		print(summary(pairscan.results, perm = pairscan.results.perm, alphas = c(0.05, 0.05, 0, 0.05, 0.05), lodcol = i, pvalues = TRUE))
		}else{
		print(summary(pairscan.results, thresholds = std.thresh, lodcol = i))	
		}	
	}
sink()

}