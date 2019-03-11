#This function adds all the derived phenotypes for 
#lean mass and fat mass. It is separated to reduce 
#clutter

add.ltm.ftm.derived.phenotypes <- function(data.obj){
	
	calc.slope <- function(pheno1.t1, pheno1.t2, pheno2.t1, pheno2.t2){
	rise <- data.obj$pheno[,pheno2.t2]-data.obj$pheno[,pheno2.t1]
	run <- data.obj$pheno[,pheno1.t2]-data.obj$pheno[,pheno1.t1]
	slope <- rise/run
	return(slope)
	}

calc.mag <- function(pheno1.t1, pheno1.t2, pheno2.t1, pheno2.t2){
	rise <- data.obj$pheno[,pheno2.t2]-data.obj$pheno[,pheno2.t1]
	run <- data.obj$pheno[,pheno1.t2]-data.obj$pheno[,pheno1.t1]
	mag <- sqrt((rise^2)+(run^2))	
	return(mag)
	}

	
	#make all the derivative phenotypes
	
	FTM1 <- data.obj$pheno[,"TTM1"] - data.obj$pheno[,"LTM1"]
	logFTM1 <- log10(FTM1)
	data.obj$pheno <- cbind(data.obj$pheno, FTM1, logFTM1)
	
	FTM2 <- data.obj$pheno[,"TTM2"] - data.obj$pheno[,"LTM2"]
	logFTM2 <- log10(FTM2)
	data.obj$pheno <- cbind(data.obj$pheno, FTM2, logFTM2)
	
	delta.ltm <- data.obj$pheno[,"LTM2"] - data.obj$pheno[,"LTM1"]
	delta.ftm <- FTM2 - FTM1
	neg.locale <- which(delta.ftm < 0)
	log.delta.ftm <- log10(abs(delta.ftm))
	log.delta.ftm[neg.locale] <- log.delta.ftm[neg.locale]*-1
	data.obj$pheno <- cbind(data.obj$pheno, delta.ltm, delta.ftm, log.delta.ftm)
	
	ltm.log.ftm.slope <- calc.slope("logFTM1", "logFTM2", "LTM1", "LTM2")
	ltm.log.ftm.mag <- calc.mag("logFTM1", "logFTM2", "LTM1", "LTM2")
	data.obj$pheno <- cbind(data.obj$pheno, ltm.log.ftm.slope, ltm.log.ftm.mag)
	
	return(data.obj)
	
}