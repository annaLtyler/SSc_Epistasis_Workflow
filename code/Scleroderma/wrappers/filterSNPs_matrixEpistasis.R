#This function uses MatrixEpistasis to filter the SSc SNPs


#code for merging cape and capeDO
#using SSc data
#make sure that if you want to run the singlescan
#there are no singlescan objects in the results directory
#otherwise, this code reads that in, and skips the singlescan
#!!! Remember to place the gene list in the results directory!!! 
testing = FALSE
exp.name <- "lung_pheno"
geno.coding <- "Dominant"

#===============================================================
# check whether we are on the cluster or home machine
#===============================================================

	cur.dir <- getwd()
	test.dir <- strsplit(cur.dir, "home")
	on.cluster = FALSE
	if(length(test.dir[[1]]) > 1){on.cluster = TRUE}
	
	if(on.cluster){		
		base.dir <- "/home/atyler/SSc"
		r.dir <- "/home/atyler/cape/"
		results.base <- "/home/atyler/"
		cape.do.dir <- "/home/atyler/code/cape"
		cape.addons.dir <- "/home/atyler/code/cape_addons"
		useful.scripts.dir <- "/home/atyler/code/useful_r_code"
		results.dir <- paste0("/home/atyler/SSc/Results/", geno.coding, "_", exp.name)
		init.dir <- "/home/atyler/SSc/data/"
		n.cores = 20
		}else{
		genotype.dir <- "~/Documents/Data/Scleroderma/GWAS_documents/geno/"
		phenotype.dir <- "~/Documents/Data/Scleroderma/GWAS_documents/pheno/phenotypes_temp"
		cape.do.dir <- "~/Documents/git_repositories/capempp"
		cape.addons.dir <- "~/Documents/git_repositories/cape_addons"
		useful.scripts.dir <- "~/Documents/git_repositories/useful_r_code"
		results.dir <- paste0("~/Documents/Data/Scleroderma/Results/", geno.coding, "_", exp.name)
		init.dir <- "~/Documents/Data/Scleroderma/project_data/"
		n.cores = 4		
		}

if(!file.exists(results.dir)){system(paste("mkdir", results.dir))}

#===============================================================
# install and load all the necessary libraries
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", "doParallel", "foreach", "caTools", "stringr", "abind", "qpcR", "biomaRt", "MatrixEpistasis")

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
all.fun <- list.files(pattern = ".R", path = cape.do.dir, full.names = TRUE)
for(i in 1:length(all.fun)){source(all.fun[i])}


all.fun <- list.files(pattern = ".R", path = cape.addons.dir, full.names = TRUE)
for(i in 1:length(all.fun)){source(all.fun[i])}


all.fun <- list.files(pattern = ".R", path = useful.scripts.dir, full.names = TRUE)
if(length(all.fun) > 0){for(i in 1:length(all.fun)){source(all.fun[i])}}

cross <- readRDS(paste0(init.dir, "cross.RData"))
geno <- readRDS(paste0(init.dir, "geno.RData"))

	
if(geno.coding == "Dominant"){
	geno[which(geno >= 0.5)] <- 1
	}
if(geno.coding == "Recessive"){
	geno[which(geno <= 0.5)] <- 0
	}

if(testing){
	sampled.idx <- sort(sample(1:length(cross$chromosome), 1000))
	cross$geno.names[[3]] <- cross$geno.names[[3]][sampled.idx]
	cross$marker.num <- cross$marker.num[sampled.idx]
	cross$chromosome <- cross$chromosome[sampled.idx]
	cross$marker.location <- cross$marker.location[sampled.idx]
	}


if(!file.exists(results.dir)){system(paste("mkdir", results.dir))}
setwd(results.dir)

data.obj <- run.cape(cross, geno, initialize.only = TRUE)


snps.to.test <- filterSNPs(data.obj, geno, pval = 0.001, chunk.size = 5000, n.cores = n.cores)
print(length(snps.to.test))


saveRDS(snps.to.test, "FilteredSNPs.RData")


