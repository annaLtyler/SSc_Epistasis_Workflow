#This is for use with binary phenotypes only
# snp1 <- "rs204990"
# snp2 <- "rs10508168"
# phenotype <- "ana"

get.snp.odds.ratio <- function(data.obj, snp1, snp2, phenotype, plot.result = FALSE){
	require("fmsb")
	
	geno <- data.obj$geno.for.pairscan
	marker.names <- unlist(lapply(strsplit(colnames(geno), "_"), function(x) x[1]))
	
	snp1.locale <- which(marker.names == snp1)
	snp2.locale <- which(marker.names == snp2)
	
	pheno.locale <- which(colnames(data.obj$pheno) == phenotype)
	pheno <- data.obj$pheno[,pheno.locale]
	

	geno.table <- cbind(geno[,snp1.locale], geno[,snp2.locale])
	geno.table <- apply(geno.table, 2, function(x) bin.vector(x, c(0, 0.5)))

	geno.pairs <- pair.matrix(c(0, 0.5), ordered = TRUE, self.pairs = TRUE)

	#===============================================================
	# internal functions
	#===============================================================
	get.combined.geno <- function(geno.table){
		combined.geno <- rep(0, nrow(geno.table))
		max.geno = max(geno.table, na.rm = TRUE)
		double.het <- apply(geno.table,1, function(x) x[1] == max.geno && x[2] == max.geno)
		combined.geno[which(double.het)] <- 1
		return(combined.geno)
		}
					
	odds.ratio <- function(exposureV, resultV){
		min.exposure = min(exposureV, na.rm = TRUE);max.exposure = max(exposureV, na.rm = TRUE)
		min.result = min(resultV, na.rm = TRUE);max.result = max(resultV, na.rm = TRUE)
		n00 <- length(intersect(which(exposureV == min.exposure), which(resultV == min.result)))
		n01 <- length(intersect(which(exposureV == min.exposure), which(resultV == max.result)))
		n10 <- length(intersect(which(exposureV == max.exposure), which(resultV == min.result)))
		n11 <- length(intersect(which(exposureV == max.exposure), which(resultV == max.result)))
		# OR <- (n00*n11)/(n01*n10)
		OR.test <- oddsratio(n11, n01, n10, n00)
		OR <- OR.test$estimate
		conf.int <- OR.test$conf.int[1:2]
		return(c(OR, conf.int))
		}
	#===============================================================

	result.mat <- matrix(NA, ncol = 3, nrow = 3)
	rownames(result.mat) <- c("OR", "min.95", "max.95")
	colnames(result.mat) <- c(snp1, snp2, "Interaction")
	
	result.mat[,1] <- odds.ratio(exposureV = geno.table[,1], resultV = pheno)
	result.mat[,2] <- odds.ratio(exposureV = geno.table[,2], resultV = pheno)
	combined.geno <- get.combined.geno(geno.table)
	result.mat[,3] <- odds.ratio(combined.geno, pheno)
	
	if(plot.result){
		par(mar = c(4,8,4,4))
		plot.new()
		plot.window(ylim = c(0.5, 3.5), xlim = c(0, max(c(as.vector(result.mat), 1), na.rm = TRUE)))
		points(y = 3:1, x = result.mat[1,], pch = 16)
		segments(y0 = 3:1, x0 = result.mat[2,], y1 = 3:1, x1 = result.mat[3,])
		par(xpd = TRUE)
		text(y = c(3:1), x = rep(-0.02,3), labels = c(snp1, snp2, "Interaction"), adj = 1)
		mtext(phenotype)
		par(xpd = FALSE)
		axis(1)
		}
	
	return(result.mat)	
	}