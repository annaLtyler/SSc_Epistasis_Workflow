#This function looks across different allele frequencies, penetrances,
#and misclassification rates to see how misclassification of patient
#status affects power to detect a SNP at different allele frequencies
#and penetrances with additive SNP coding

snp.power <- function(num.pt = 833, num.exp = 100, sig.level = 0.05){
	
	maf <- seq(0.05, 0.5, 0.05)
	# pen <- seq(1, 100, 1)
	pen <- 0.25
	mis.class <- seq(0,25,1)
	
	# p.array <- array(NA, dim = c(length(maf), length(pen), length(mis.class)))
	# power.array <- array(0, dim = c(length(maf), length(pen), length(mis.class)))

run.exp <- function(e){

	p.array <- array(NA, dim = c(length(maf), length(pen), length(mis.class)))
	dimnames(p.array) <- list(maf, pen, mis.class)
	# cat("\n\nExperiment", e, "\n")
	for(m in 1:length(maf)){
		# cat("\nMinor Allele Frequency:", maf[m], "\n")
		#simulate a SNP with a given allele frequency
		snp <- matrix(0, nrow = num.pt, ncol = 2)
		num.alleles <- num.pt*2
		num.ma <- round(num.alleles*maf[m])
		allele.placement <- sample(1:num.alleles, num.ma)
		snp[allele.placement] <- 1
		geno <- apply(snp, 1, mean)
		
		#calculate case control status based on penetrance
		#each allele contributed an additive penetrance component
		for(p in 1:length(pen)){
			# report.progress(p, length(pen))
			status <- rep(0, num.pt)
			risk <- geno*pen[p]
			prob <- runif(num.pt, 0, 1)
			with.disease <- which(apply(cbind(risk,prob), 1, function(x) x[1] > x[2]))
			if(length(with.disease) > 0){
				status[with.disease] <- 1
				}
			
			#now flip some bits to represent misclassification
			for(mc in 1:length(mis.class)){
				# cat("\tMis-classification:", pen[mc], "\n")
				num.mis <- round((num.pt*mis.class[mc])/100)
				if(num.mis > 0){
					wrong.status <- sample(1:num.pt, num.mis)
					status[which(status[wrong.status] == 1)] <- 0
					status[which(status[wrong.status] == 0)] <- 1
					}		
			#run the model to see how significant the SNP is
			model <- lm(status~geno)
			p.array[m,p,mc] <- summary(model)$coefficients[2,4]
			} #end looping through mis classification rates
		} #end looping through penetrances
	}# end looping through minor allele frequencies

	# sig.locale <- which(p.array <= sig.level)
	# power.array[sig.locale] <- power.array[sig.locale] + 1

	return(p.array)
} #end run experiment function

		sfInit(parallel = TRUE, cpus = 4)
		sfExport("p.array", "run.exp", "maf", "pen", "mis.class", "num.pt", "num.exp", "sig.level")
		p.arrays <- sfLapply(1:num.exp, run.exp)
		sfRemove("p.array", "run.exp", "maf", "pen", "mis.class", "num.pt", "num.exp", "sig.level")
		sfStop()


	saveRDS(p.array, "pValueArray.RData")
	saveRDS(power.array, "powerArray.RData")

	#plot p values for different parameters
	plot.array <- function(data.array, type = c("pval", "power")){
		val.mat <- array(NA, dim = c(dim(data.array[[1]])[1], dim(data.array[[1]])[3], num.exp))
		for(e in 1:length(data.array)){
			if(type == "pval"){
				exp.array <- log10(data.array[[e]])*-10
				ylab = "-log10 p values"
				}
			if(type == "power"){
				exp.array <- data.array[[e]]/num.exp
				ylab = "Power"
				}
			
			val.mat[,,e] <- exp.array[,1,]
			}
		
		plot.new()
		plot.window(xlim = c(0,26), ylim = c(min(val.mat), max(val.mat)))
		test <- apply(val.mat, 3, function(x) boxplot(x, add = TRUE, axes = FALSE))
		axis(1); axis(2)

		plot.new()
		plot.window(xlim = c(0,10), ylim = c(min(val.mat), max(val.mat)))
		test <- apply(val.mat, 3, function(x) boxplot(t(x), add = TRUE, axes = FALSE))
		axis(1); axis(2)

		}

		

	pdf("SNP.Power.Values.pdf")
	plot.array(p.array, 1, maf, type = "pval")
	plot.array(p.array, 2, pen, type = "pval")
	plot.array(p.array, 3, mis.class, type = "pval")
	dev.off()
	

	pdf("SNP.Power.Values.pdf")
	plot.array(power.array, 1, maf, type = "power")
	plot.array(power.array, 2, pen, type = "power")
	plot.array(power.array, 3, mis.class, type = "power")
	dev.off()
	



}