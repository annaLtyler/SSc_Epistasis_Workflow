#code for merging cape and capeDO
#using SSc data
#make sure that if you want to run the singlescan
#there are no singlescan objects in the results directory
#otherwise, this code reads that in, and skips the singlescan

test.genes <- c("FHIT", "GPSM3", "NOTCH4")
upstream.buffer = downstream.buffer = 1000
geno.coding = "Dominant"
#===============================================================
# check whether we are on the cluster or home machine
#===============================================================

	cur.dir <- getwd()
	test.dir <- strsplit(cur.dir, "home")
	on.cluster = FALSE
	if(length(test.dir[[1]]) > 1){on.cluster = TRUE}
	
	if(on.cluster){		
		base.dir <- "/home/atyler/SSc"
		cape.do.dir <- "/home/atyler/code/cape/"
		cape.addons.dir <- "/home/atyler/code/cape_addons"
		useful.scripts.dir <- "/home/atyler/code/useful_r_code/"
		results.dir <- "/home/atyler/SSc/Results/Dominant_max.pair_4ET_fine"
		n.cores = 20
		}else{
		cape.do.dir <- "~/Documents/git_repositories/capempp"
		cape.addons.dir <- "~/Documents/git_repositories/cape_addons"
		useful.scripts.dir <- "~/Documents/git_repositories/useful_r_code/"
		data.dir <- "~/Documents/Data/Scleroderma/Results/Dominant_max.pair_4ET"
		results.dir <- "~/Documents/Data/Scleroderma/Results/Dominant_max.pair_4ET_fine"
		geno.dir <- "~/Documents/Data/Scleroderma/GWAS_documents/geno"
		n.cores = 4
		}

if(!file.exists(results.dir)){system(paste("mkdir", results.dir))}

#===============================================================
# install and load all the necessary libraries
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", "doParallel", "foreach", "caTools", "stringr", "abind", "qpcR")

if(on.cluster){
	all.packages <- installed.packages(lib.loc = "/home/atyler/R/x86_64-pc-linux-gnu-library/3.4")
	
	package.locale <- match(needed.packages, all.packages[,1])
	missing.packages <- needed.packages[which(is.na(package.locale))]
	
	if(length(missing.packages) > 0){install.packages(missing.packages, lib = "/home/atyler/R/x86_64-pc-linux-gnu-library/3.4", repos = "https://cloud.r-project.org")}
	}

for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}
#===============================================================


#===============================================================
#source all the code from needed packages
#===============================================================
setwd(cape.do.dir)
all.fun <- list.files(pattern = ".R")
for(i in 1:length(all.fun)){source(all.fun[i])}

setwd(cape.addons.dir)
all.fun <- list.files(pattern = ".R")
for(i in 1:length(all.fun)){source(all.fun[i])}

setwd(useful.scripts.dir)
all.fun <- list.files(pattern = ".R")
if(length(all.fun) > 0){for(i in 1:length(all.fun)){source(all.fun[i])}}

#===============================================================
if(!on.cluster){
	library(biomaRt)
	hum <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	# hum <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")


	gene.location <- getBM(c("external_gene_name", "chromosome_name", "start_position", "end_position"), "external_gene_name", test.genes, hum)
	#This causes a NAs introduced by coercion warning, but that's fine
	chr.loc <- which(!is.na(as.numeric(gene.location[,"chromosome_name"])))
	gene.location <- gene.location[chr.loc,]


	setwd(data.dir)
	cross <- readRDS("cross.RData")
	geno <- readRDS(paste0(geno.dir, "/geno.RData"))
	snp.gene.info <- as.matrix(read.table("~/Documents/Data/Scleroderma/GWAS_documents/geno/SNP.Gene.Info.txt", sep = "\t", stringsAsFactors = F))
	gene.snp.list <- vector(mode = "list", length = length(test.genes))
	names(gene.snp.list) <- test.genes
	for(g in 1:length(test.genes)){
		gene.name <- test.genes[g]
		gene.locale <- which(snp.gene.info[,2] == test.genes[g])
		gene.snps <- snp.gene.info[gene.locale,1]
		snps.to.test <- intersect(cross$geno.names[[3]], gene.snps)
		gene.snp.list[[g]] <- as.vector(snps.to.test)
		}
	
	snp.list <- unlist(gene.snp.list)
	
	if(on.cluster){
		geno <- geno$geno
		}
	
	ind <- match(rownames(cross$pheno), rownames(geno))
	new.geno.for.pairscan <- geno[ind,2,snp.list]
		
	if(geno.coding == "Dominant"){
		new.geno.for.pairscan[which(new.geno.for.pairscan >= 0.5)] <- 1
		}
	if(geno.coding == "Recessive"){
		new.geno.for.pairscan[which(geno <= new.geno.for.pairscan)] <- 0
		}
	

	#===============================================================
	#keep track of the experimental parameters
	#===============================================================
	setwd(results.dir)
	
	cross$geno.for.pairscan <- new.geno.for.pairscan
		
	cross.pairscan <- pairscan(data.obj = cross, geno.obj = geno, scan.what = "ET", total.perm = 0, min.per.geno = NULL, max.pair.cor = 0.5, verbose = TRUE, num.pairs.limit = Inf, overwrite.alert = FALSE, run.parallel = TRUE, n.cores = n.cores)
		
	saveRDS(cross.pairscan, "cross.pairscan.RData")
	saveRDS(cross, "cross.RData")
			
	orig.pairscan <- readRDS(paste0(data.dir, "/cross.pairscan.RData"))
	
	cross.pairscan$pairscan.perm <- orig.pairscan$pairscan.perm
	cross.pairscan$pairs.tested.perm <- orig.pairscan$pairs.tested.perm
	saveRDS(cross.pairscan, "cross.pairscan.RData")
	}else{
	setwd(results.dir)
	cross <- readRDS("cross.RData")
	cross.pairscan <- readRDS("cross.pairscan.RData")	
	}

#===============================================================
# run reprametrization
#===============================================================
	
cat("Performing error propagation on cape coefficients...\n")
cross <- error.prop(data.obj = cross, pairscan.obj = cross.pairscan, perm = FALSE, verbose = TRUE, n.cores = n.cores)
	
saveRDS(cross, "cross.RData")
	
cat("Calculating p values...\n")
cross <- calc.p(data.obj = cross, pval.correction = "fdr", verbose = TRUE, run.parallel = TRUE, n.cores = n.cores)
	
	cat("Calculating directed influences...\n")
	cross <- direct.influence(data.obj = cross, pairscan.obj = cross.pairscan, transform.to.phenospace = TRUE, verbose = TRUE, pval.correction = "fdr", save.permutations = TRUE, n.cores = n.cores)
	
	saveRDS(cross, "cross.RData")
	
	pdf("variant.influences.pdf", width = 10, height = 7)
	inf.mat <- plotVariantInfluences(data.obj = cross, p.or.q = 0.05, standardize = FALSE, not.tested.col = "lightgray", covar.width = 30, pheno.width = 30)
	dev.off()
	
	writeVariantInfluences(data.obj = cross, p.or.q = 0.2, filename = "Variant.Influences.csv")
	
	# # cat("Generating the network...\n")
	cross <- get.network(cross, p.or.q = 0.05, collapse.linked.markers = FALSE)
	cross <- get.network(cross, p.or.q = 0.05, threshold.power = 1, collapse.linked.markers = TRUE, plot.linkage.blocks = TRUE)
	
	saveRDS(cross, "cross.RData")
	
	if(!on.cluster){
	write.data.for.shiny(cross, "SSc_FHIT_NOTCH4_GPSM3", path = "~/Documents/R/shiny/capeDO_interactions/data")
	}