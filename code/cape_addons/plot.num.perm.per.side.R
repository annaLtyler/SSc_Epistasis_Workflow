

plot.num.perm.per.side <- function(data.obj, min.perm = 1000, max.perm = "max", num.bins = 100){

	cur.dir <- getwd()
	r.dir <- "~/Documents/git_repositories/cape"
	setwd(r.dir)
	all.fun <- list.files(pattern = "*.R")
	for(i in 1:length(all.fun)){source(all.fun[i])}
	setwd(cur.dir)


	influences.org <- data.obj$var.to.var.influences
	influences.perm <- data.obj$var.to.var.influences.perm

	
	mat12.perm <- as.numeric(influences.perm[,3]) / as.numeric(influences.perm[,4])
	mat21.perm <- as.numeric(influences.perm[,5]) / as.numeric(influences.perm[,6])
	mat12.mat21.perm <- c(mat12.perm, mat21.perm)

	if(max.perm == "max"){
		max.perm <- length(mat12.mat21.perm)
		}
	perm.seq <- round(segment.region(min.perm, max.perm, num.bins, "ends"))


	num.neg <- rep(NA, length(perm.seq))
	num.pos <- rep(NA, length(perm.seq))
	for(i in 1:length(perm.seq)){
		chopped.null <- mat12.mat21.perm[1:perm.seq[i]]
		num.neg[i] <- length(which(chopped.null <= median(chopped.null)))
		num.pos[i] <- length(which(chopped.null > median(chopped.null)))
		}

	pdf("pval.resolution.pdf")	
	max.y <- max(c(log10(1/num.neg)*-1, log10(1/num.pos)*-1, log10(1/perm.seq)*-1))
	min.y <- min(c(log10(1/num.neg)*-1, log10(1/num.pos)*-1, log10(1/perm.seq)*-1))
	plot(perm.seq, log10(1/num.neg)*-1, col = "blue", lwd = 3, type = "l", ylim = c(min.y, max.y), xlab = "Total Permutations", ylab = "-log10 minimum p value")
	points(perm.seq, log10(1/perm.seq)*-1, col = "black", lwd = 3, type = "l")
	legend("bottomright", legend = c(paste("total permutations:", signif(log10(1/max(perm.seq))*-1, 2)), paste("half permutations:", signif(log10(1/max(num.neg))*-1, 2))), col = c("black", "blue"), lty = 1, lwd = 3)
	dev.off()

	}