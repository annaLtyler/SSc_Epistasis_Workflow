#This function takes cape data and geno objects
#and writes out the files required for R/qtl2
#input. This is for reading in the data as an
#R/qtl2 object for calculating kinship matrices

write_cape_as_qtl2 <- function(data.obj, geno.obj, scan.what = c("eigentraits", "raw.traits"), min.geno.prob = 0.4, cross.type = "do", filename = "project.zip", project.dir = ".", map.type = c("physical", "genetic"), convert.bp.to.Mbp = TRUE){
	
	map.type = map.type[1]
	cross.type = cross.type[1]
	project.name <- strsplit(filename, ".zip")[[1]][1]
	
	cat("Writing phenotype data...\n")
	pheno <- get.pheno(data.obj, scan.what)
	pheno.table <- cbind(rownames(pheno), pheno)
	colnames(pheno.table)[1] <- "ID"
	write.table(pheno.table, "project_pheno.csv", sep = ",", quote = FALSE, row.names = FALSE)
	control.line.pheno <- "pheno: project_pheno.csv"
	
	cat("Writing covariate data...\n")
	covar.info <- get.covar(data.obj)
	covar.table <- cbind(rownames(covar.info$covar.table), covar.info$covar.table)
	colnames(covar.table) <- c("ID", covar.info$covar.names)
	write.table(covar.table, "project_covar.csv", sep = ",", quote = FALSE, row.names = FALSE)
	control.line.covar <- "covar: project_covar.csv"
	
	cat("Converting genotype probabilities to haplotypes...\n")
	geno.table <- genoprob2genotype(geno.obj, min.geno.prob)
	geno.table <- cbind(rownames(geno.table), geno.table)
	colnames(geno.table)[1] <- "ID"
	write.table(geno.table, "project_geno.csv", sep = ",", quote = FALSE, row.names = FALSE)
	control.line.geno <- "geno: project_geno.csv"
	
	cat("Writing map data...\n")
	map.table <- cbind(data.obj$geno.names[[3]], data.obj$chromosome, data.obj$marker.location)
	colnames(map.table) <- c("marker", "chr", "pos")
	if(convert.bp.to.Mbp){
		map.table[,3] <- as.numeric(map.table[,3])/1e6
		}
	write.table(map.table,  "project_pmap.csv", sep = ",", quote = FALSE, row.names = FALSE)	
	
	if(map.type == "physical"){
		control.line.map <- "pmap: project_pmap.csv"
		}else{
		control.line.map <- "gmap: project_pmap.csv"	
		}

	alleles <- colnames(geno.obj)
	allele.table <- cbind(rep(" ", length(alleles)), rep("-", length(alleles)), alleles)

	genotypes <- pair.matrix(alleles, self.pairs = TRUE)
	genotype.table <- apply(cbind(genotypes, 1:nrow(genotypes)), 1, function(x) paste(x[1], x[2], ": ", x[3], collapse = "", sep = ""))

	cat("Writing control file and compiling project...\n")
	write.table("# project description", "project.yaml", quote = FALSE,row.names = FALSE, col.names = FALSE)
	write.table(paste("crosstype:", cross.type), "project.yaml", append = TRUE, quote = FALSE,row.names = FALSE, col.names = FALSE)
	write.table(control.line.geno, "project.yaml", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(control.line.pheno, "project.yaml", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(control.line.covar, "project.yaml", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(control.line.map, "project.yaml", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table("alleles:", "project.yaml", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(allele.table, "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table("genotypes:", "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(genotype.table, "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table("sex:", "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(" covar: Sex", "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(" 0: female", "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table(" 1: male", "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table("na.strings:", "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table("- NA", "project.yaml", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	
	project.files <- list.files(pattern = "project")
	system(paste("mkdir", project.name))
	for(i in 1:length(project.files)){
		system(paste("mv", project.files[i], project.name))
		}
	system(paste("zip -r", filename, project.name))
	system(paste("rm -r", project.name))

}