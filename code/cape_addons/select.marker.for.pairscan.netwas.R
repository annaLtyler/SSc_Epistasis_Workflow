#This script selects markers for the pair scan.
#to specify individual markers, create a list
#called specific markers. The names of the list
#elements should be the loci and each element
#should be the alleles wanted for that marker
#non-allelic covar should specify which covariates
#do not have alleles. These will be pared down
#so that there is only one representative. 
#if num.markers is specified, a threshold is 
#set such that the given number of markers 
#(within the tolerance specifiec) are selected.
#testing: snp.info.file = "SNP.Gene.Info.txt"; netwas.sig.val = 0.01; tissue.type = "skin"; top.netwas.genes = 100; select.by = "mean"; plot.diagnostics = TRUE; alpha.thresh = NULL; t.thresh = NULL; specific.markers = NULL; num.markers = 100; step.size = NULL; tolerance = 10; verbose = TRUE

select.markers.for.pairscan.netwas <- function(data.obj, geno.obj = NULL, singlescan.obj, snp.info.file = "SNP.Gene.Info.txt", netwas.sig.val, tissue.type, top.netwas.genes = NULL, select.by = c("mean", "median", "max"), plot.diagnostics = TRUE, alpha.thresh = NULL, t.thresh = NULL, specific.markers = NULL, num.markers = NULL, step.size = NULL, tolerance = 10, verbose = TRUE){
	
	
	data.obj1 <- select.markers.for.pairscan(data.obj, geno.obj = geno.obj, singlescan.obj, alpha.thresh = alpha.thresh, t.thresh = t.thresh, specific.markers = specific.markers, num.markers = num.markers, step.size = step.size, tolerance = tolerance, verbose = FALSE)

	if(verbose){cat(dim(data.obj1$geno.for.pairscan)[2], "SNPs have been selected based on effect size.\n")}

	data.obj2 <- filter.snps.netwas(data.obj, geno.obj = geno.obj, snp.info.file = snp.info.file, sig.val = netwas.sig.val, tissue.type = tissue.type, top.n.genes = top.netwas.genes, select.by = select.by, verbose = FALSE, plot.diagnostics = plot.diagnostics)

	if(verbose){cat(dim(data.obj2$geno.for.pairscan)[2], "SNPs have been selected based NetWAS results.\n")}

	
	geno.mat1 <- data.obj1$geno.for.pairscan
	geno.mat2 <- data.obj2$geno.for.pairscan
	

	snps1 <- colnames(geno.mat1)
	snps2 <- colnames(geno.mat2)
	
	full.snp.list <- unique(c(snps1, snps2))
	
	if(verbose){cat("There are", length(full.snp.list), "unique SNPs between the two lists.\n")}
	
	snp.num <- get.marker.num(data.obj, full.snp.list)	
	snp.order <- order(snp.num)
	
	ordered.snps <- full.snp.list[snp.order]
	
	geno <- get.geno(data.obj, geno.obj)
	
	comb.geno.for.pairscan <- geno[,ordered.snps]
	lin.ind.geno <- get.linearly.independent(geno.matrix = comb.geno.for.pairscan)
	
	if(verbose){cat(length(lin.ind.geno$rejected.markers), "SNPs were rejected for being non-linearly independent.\n")}
	
	data.obj$geno.for.pairscan <- lin.ind.geno$independent.markers
	
	if(verbose){cat("There are", dim(lin.ind.geno$independent.markers)[2], "SNPs in the genotype matrix for the pairscan.\n")}
	
	return(data.obj)
	
	
}