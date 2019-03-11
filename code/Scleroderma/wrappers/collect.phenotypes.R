labs <- as.matrix(read.table("~/Documents/Grants/R21/Scleroderma/GWAS_documents/phenotypes_temp/phs000357.v1.pht002331.v1.p1.c2.Systemic_Sclerosis_Lab_Results.SSCAA.txt", sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE))
phenotypes <- c("dbGaP.SampID", "SAMPID", "ana", "centromere", "nucleolar", "speckled", "topo", "rna_ind")
lab.codes <- c(NA, NA, rep("pnn", 6))
phenotype.locale <- match(phenotypes, colnames(labs))
lab.table <- labs[,phenotype.locale]

clinical <- as.matrix(read.table("~/Documents/Grants/R21/Scleroderma/GWAS_documents/phenotypes_temp/phs000357.v1.pht002332.v1.p1.c2.Systemic_Sclerosis_Clinical_Attributes.SSCAA.txt", sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE))
phenotypes <- c("subjid", "lv_crei", "hr_srci")
clinical.codes <- c(NA, "ann", "yn")
phenotype.locale <- match(phenotypes, colnames(clinical))
clinical.table <- clinical[,phenotype.locale]

demo <- as.matrix(read.table("~/Documents/Grants/R21/Scleroderma/GWAS_documents/phenotypes_temp/phs000357.v1.pht002333.v1.p1.c2.Systemic_Sclerosis_Demographics.SSCAA.txt", sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE))
phenotypes = c("subjid", "gender", "ssctype", "dis_diagag", "rayonset", "nronsetage")
demo.codes <- c(NA, "sex", NA, rep("num", 3))
phenotype.locale <- match(phenotypes, colnames(demo))
demo.table <- demo[,phenotype.locale]

lung <- as.matrix(read.table("~/Documents/Grants/R21/Scleroderma/GWAS_documents/phenotypes_temp/phs000357.v1.pht002334.v1.p1.c2.Systemic_Sclerosis_Pulmonary_Attributes.SSCAA.txt", sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE))
phenotypes <- c("subjid", "hp_ildi", "hc_pahi", "ec_phi", "ec_phed", "ec_lvei", "ec_pei", "pf1_date", "pf1_fvcp", "pf1_dlcp", "pf1_fevp", "pf1_tlcp", "pf2_date", "pf2_fvcp", "pf2_dlcp", "pf2_fevp", "pf2_tlcp")
lung.codes <- c(NA, rep("yn", 3), "num", rep("yn", 2), rep("num", 10))
phenotype.locale <- match(phenotypes, colnames(lung))
lung.table <- lung[,phenotype.locale]

#changes positive/negative/null to 1/0/NA
pnn.code <- function(entryV){
	entryV[which(entryV == 1)] <- 1
	entryV[which(entryV == 2)] <- 0
	entryV[which(entryV == 9)] <- NA
	return(entryV)
	}

yn.code <- function(entryV){
	entryV[which(entryV == "N")] <- 0
	entryV[which(entryV == "Y")] <- 1
	entryV[which(entryV == "Z")] <- NA
	return(entryV)
	}

#changes abnormal/normal/null to 1/0/NA
ann.code <- function(entryV){
	entryV[which(entryV == "A")] <- 1
	entryV[which(entryV == "N")] <- 0
	entryV[which(entryV == "Z")] <- NA
	return(entryV)
	}
	
sex.code <- function(entryV){
	entryV[which(entryV == "Male")] <- 1
	entryV[which(entryV == "Female")] <- 0
	entryV[which(entryV == "Unknown")] <- NA
	return(entryV)
	}
	
num.code <- function(entryV){
	entryV[which(entryV == "9999")] <- NA
	entryV[which(entryV == "999.9")] <- NA
	return(as.numeric(entryV))
	}
	
recode.ssc <- function(entryV){
	new.table <- matrix(0, nrow = length(entryV), ncol = 3)
	colnames(new.table) <- c("limited", "diffuse", "sine")
	new.table[which(entryV == 1),1] <- 1
	new.table[which(entryV == 2),2] <- 1
	new.table[which(entryV == 6),3] <- 1
	return(new.table)
	}
	
recode.table <- function(orig.table,codes){
	coded.table <- orig.table
	for(i in 1:dim(orig.table)[2]){
	if(!is.na(codes[i])){
		if(codes[i] == "ann"){
			coded.table[,i] <- ann.code(orig.table[,i])
			}
		if(codes[i] == "pnn"){
			coded.table[,i] <- pnn.code(orig.table[,i])
			}
		if(codes[i] == "sex"){
			coded.table[,i] <- sex.code(orig.table[,i])
			}
		if(codes[i] == "num"){
			coded.table[,i] <- num.code(orig.table[,i])
			}
		if(codes[i] == "yn"){
			coded.table[,i] <- yn.code(orig.table[,i])
			}

			}
		}
	return(coded.table)
	}
	
#recode all values (except SSC type)
coded.lab <- apply(recode.table(lab.table, lab.codes), 2, as.numeric)
coded.demo <- apply(recode.table(demo.table, demo.codes), 2, as.numeric)
coded.clinical <- apply(recode.table(clinical.table, clinical.codes), 2, as.numeric)
coded.lung <- apply(recode.table(lung.table, lung.codes), 2, as.numeric)

#recode SSc type
ssc.type <- recode.ssc(coded.demo[,"ssctype"])
coded.demo <- cbind(coded.demo[,-3], ssc.type)

#put the tables together making sure to match up patient numbers
all.pt <- unique(c(coded.lab[,"SAMPID"], coded.demo[,"subjid"], coded.clinical[,"subjid"], coded.lung[,"subjid"]))


pt.locale <- match(all.pt, coded.lab[,"SAMPID"])
final.table <- coded.lab[pt.locale,]
pt.locale <- match(all.pt, coded.demo[,"subjid"])
final.table <- cbind(final.table, coded.demo[pt.locale,])
pt.locale <- match(all.pt, coded.clinical[,"subjid"])
final.table <- cbind(final.table, coded.clinical[pt.locale,])
pt.locale <- match(all.pt, coded.lung[,"subjid"])
final.table <- cbind(final.table, coded.lung[pt.locale,])
id.col <- which(colnames(final.table) == "subjid")
final.table <- final.table[,-id.col]

write.table(final.table, "SSc_phenotypes.csv", quote = FALSE, sep = ",", row.names = FALSE)
num.pts <- apply(final.table, 2, function(x) length(which(!is.na(x))))
write.table(num.pts, "SSC_phenotype_numbers.csv", sep = ",", quote = FALSE, col.names = FALSE)

