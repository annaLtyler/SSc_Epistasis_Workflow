#This function calculates risk ratios for pairs of SNPs
#as well as a number of interaction terms as defined in
#Kalilani and Atashili 2006
#the outcome threshold is for when we have corrected the
#phenotype for a covariate and the outcome is no longer
#0 or 1 

get.snp.int.risk.ratio <- function(data.obj, geno.obj = NULL, snp1, snp2, phenotype, allele.which = 2, covar = NULL, trait.type = c("ET", "Normalized", "Raw"), outcome.threshold = 0.5, geno.coding = c("Additive", "Dominant", "Recessive"), plot.results = TRUE, las = 2){
	
	if(is.null(geno.obj)){
		geno <- data.obj$geno.for.pairscan
		snp.locale <- match(c(snp1, snp2), colnames(geno))
		geno.table <- geno[,snp.locale]
		}else{
		geno <- get.geno(data.obj, geno.obj)
		geno.table <- geno[,allele.which,c(snp1, snp2)]
		}

	cor.pheno <- get.pheno(data.obj, trait.type[1], covar)
	pheno.locale <- which(colnames(data.obj$pheno) == phenotype)
	pheno <- cor.pheno[,pheno.locale]
	
	geno.coding = geno.coding[1]
	if(geno.coding == "Dominant"){
		geno.table[which(geno.table >= 0.5)] <- 1
		}
	if(geno.coding == "Recessive"){
		geno.table[which(geno.table <= 0.5)] <- 0
		}					
	
	#calculate risk ratio for each SNP
	snp1.rr <- risk.ratio(exposure.vector = geno.table[,1], outcome.vector = pheno, outcome.threshold = outcome.threshold)
	snp2.rr <- risk.ratio(exposure.vector = geno.table[,2], outcome.vector = pheno, outcome.threshold = outcome.threshold)
	
	#calculate the interaction contrast (IC)
	risks <- risk.interactions(exposure.vector1 = geno.table[,1], exposure.vector2 = geno.table[,2], outcome.vector = pheno, plot.results = FALSE, outcome.threshold = outcome.threshold)
	IC <- as.numeric(risks[5])

	R11 <- risks["1/1"]/risks["0/0"]
	R01 <- risks["0/1"]/risks["0/0"]
	R10 <- risks["1/0"]/risks["0/0"]
	
	# c(R01, R10, R11)
	if(plot.results){
		barplot(c(R10, R01, R11), names = c(snp1, snp2, "both"), ylab = "Risk Relative to Neither SNP", main = phenotype, ylim = c(0, max(c(R01, R10, R11, 1))), las = las)
		abline(h = 1)
		}

	#calculate different interaction terms for double exposures
	mult.int <- as.numeric(R11/(R01*R10))
	
	#ICR: "the excess risk due to interaction 
	#relative to the risk without exposure"
	#There is no interaction effect if the value is 0
	ICR <- as.numeric(IC/risks["0/0"])
	# ICR.check <- (risks["1/1"]/risks["0/0"]) - (risks["0/1"]/risks["0/0"]) - (risks["1/0"]/risks["0/0"]) + 1
	# ICR.check <- R11 - R01 - R10 + 1
	
	#The relative excess risk due to interaction
	RERI <- as.numeric(R11 - R10 - R01 + 1)
	
	#amount total effect attributable to SNP1 alone
	snp1.effect <- as.numeric((R10 - 1)/(R11 - 1))

	#amount of total effect attributable to SNP2 alone
	snp2.effect <- as.numeric((R01 - 1)/(R11 - 1))
	
	#amount of total effect attributable to the interaction
	int.effect <- as.numeric(RERI/(R11 - 1))
	
	#sum(c(snp1.effect, snp2.effect, int.effect))
	
	#AP: "the attributable proportion of disease 
	#which is due to interaction among persons 
	#with both exposures"
	#There is no interaction effect if the value is 0
	AP <- as.numeric(IC/risks["1/1"])
	# AP.check <- ICR/(risks["1/1"]/risks["0/0"])
	# AP.check <- ICR/R11

	#S: "the excess risk from double exposure 
	#with an interaction, relative to the 
	#excess risk from double exposure 
	#when there is no interaction."
	#There is no interaction effect if the value is 1
	#This measure can be difficult to interpret if one allele is
	#causative and the other is preventative (VanderWeele and Knol 2014)
	S <- as.numeric((risks["1/1"] - risks["0/0"])/(risks["1/0"] - risks["0/0"] + risks["0/1"] - risks["0/0"]))
	# S.check <- (R11 - 1)/(R10+R01-2)
	# S.check <- (R11 - 1)/(R10 - 1 + R01 - 1)
	
	results <- matrix(c(snp1.rr, 
						snp2.rr, 
						R10, 
						R01, 
						R11, 
						IC, 
						mult.int, 
						RERI, 
						AP, 
						S), 
						ncol = 1)
						
	rownames(results) <- c(paste0(snp1, ".risk ratio"), 
						paste0(snp2, ".risk ratio"), 
						paste0("RR.", snp1, ".to.no.SNPs"), 
						paste0("RR.", snp2, ".to.no.SNPs"),
						"RR.both.SNPs.to.no.SNPs",
						"add.IC", 
						"multi.int", 
						"RERI.or.ICR", 
						"AP", 
						"S")

	return(results)	
	}