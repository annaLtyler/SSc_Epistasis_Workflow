build.SSc.cape.data <- function(){
	
	pheno <- as.matrix(read.csv("~/Documents/Data/Scleroderma/GWAS_documents/pheno/phenotypes_temp/SSc_phenotypes.csv", stringsAsFactors = FALSE))

	demo.info <- as.matrix(read.table("~/Documents/Data/Scleroderma/GWAS_documents/pheno/phenotypes_temp/phs000357.v1.pht002333.v1.p1.c2.Systemic_Sclerosis_Demographics.SSCAA.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
	sex <- demo.info[,"gender"]
	sex[which(str_trim(sex) == "Female")] <- 0
	sex[which(str_trim(sex) == "Male")] <- 1

	new.pheno <- cbind(pheno, sex)
	new.pheno <- apply(new.pheno, 2, as.numeric)
	rownames(new.pheno) <- pheno[,2]

	head(new.pheno)
	
	geno <- readRDS('~/Documents/Data/Scleroderma/GWAS_documents/geno/GenotypeMatrixRecoded.RData')
	rownames(geno) <- rownames(new.pheno)

	snp.info <- as.matrix(read.table("~/Documents/Data/Scleroderma/GWAS_documents/geno/dbGaPcases_newID2.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE))

	# identical(snp.info[,2], colnames(geno))
	
	#take out the markers on the sex chromosomes
	x.locale <- which(snp.info[,1] == 23)
	snp.info <- snp.info[-x.locale,]
	geno <- geno[,-x.locale]
	
	#filter by minor allele frequency
	maf.thresh <- 0.1
	all.maf <- unlist(lapply_pb(1:ncol(geno), function(x) mean(geno[,x], na.rm = TRUE)))
	#hist(all.maf)
	to.keep <- which(all.maf >= maf.thresh)
	geno <- geno[,to.keep]
	snp.info <- snp.info[to.keep,]
	
	geno.obj <- list()
	geno.obj$geno <- geno
	geno.obj$geno.names <- dimnames(geno)

	#===============================================================
	# build the cape object
	#===============================================================
	cross <- list()
	cross$pheno <- new.pheno
	cross$chromosome <- as.numeric(snp.info[,1])
	cross$marker.names <- snp.info[,2]
	cross$marker.location <- as.numeric(snp.info[,4])
	cross$marker.num <- 1:nrow(snp.info)
	
	cross.obj <- cape2mpp(cross, geno.obj)
	cross <- cross.obj$data.obj
	geno.obj <- cross.obj$geno.obj

	
#===============================================================
	
	cross.geno <- delete.underscore(data.obj = cross, geno.obj = geno.obj)
	cross <- cross.geno$data.obj
	geno <- cross.geno$geno.obj
	cross <- remove.unused.markers(cross, geno)

	#===============================================================
	saveRDS(cross, "~/Documents/Data/Scleroderma/project_data/cross.init.test.RData")
	saveRDS(geno, "~/Documents/Data/Scleroderma/project_data/geno.test.RData")

	}