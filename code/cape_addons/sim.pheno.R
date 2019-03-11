#This function simulates phenotypes that are a combination of linearly separable and not linearly separable
#You can add QTL by adding a genotype object and specifying the number of QTL to be generated.

sim.pheno <- function(num.pheno = 5, num.lin.sep = 3, num.ind = 100, cont.per.non.ind = 2, rel.cont = c(0.5, 0.5), geno.obj = NULL, num.qtl = 0){
	
	if(num.qtl > 0 && is.null(geno.obj)){
		stop("A genotype matrix must be supplied to generate QTL.")
		}

	num.non.ind <- num.pheno - num.lin.sep
	
	ind.pheno.mat <- matrix(rnorm(num.lin.sep*num.ind), ncol = num.lin.sep, nrow = num.ind)
	
	if(num.non.ind > 0){
		#=================================================================	
		#generate a random matrix containing the contributions of each
		#independent phenotype to each non-independent phenotype
		#=================================================================	
		cont.mat <- matrix(0, nrow = num.non.ind, ncol = num.lin.sep)
		samples.used <- matrix(NA, nrow = num.non.ind, ncol = cont.per.non.ind)
	
		for(i in 1:num.non.ind){
			cont.which <- sample(1:num.lin.sep, cont.per.non.ind)
			samples.used[i,] <- cont.which
			if(i > 1){
				used.before <- 1
				while(sum(used.before) > 0){
					used.before <- apply(samples.used[1:(i-1),,drop=FALSE], 1, function(x) identical(x, cont.which))
					if(any(used.before)){
						used.before <- 1
						cont.which <- sample(1:num.lin.sep, cont.per.non.ind)
						samples.used[i,] <- cont.which
						}
					}
				}
			cont.mat[i,cont.which] <- sample(rel.cont)
			}
		#=================================================================	
		
		#=================================================================	
		#generate the non-independent phenotypes
		#=================================================================	
		non.ind.pheno.mat <- matrix(NA, ncol = num.non.ind, nrow = num.ind)
		for(i in 1:num.non.ind){
			#This phenotype will be determined by the phenotypes given
			#in cont.mat
			non.ind.pheno.mat[,i] <- cont.mat[i,] %*% t(ind.pheno.mat)		
			}
	
		full.pheno <- cbind(ind.pheno.mat, non.ind.pheno.mat)
		
		#add a little noise to all the phenotypes
		noise = rnorm((nrow(full.pheno)*ncol(full.pheno)), mean = 0, sd = 0.05)
		full.pheno <- full.pheno + noise
		}else{
			full.pheno  <- ind.pheno.mat
			}
	colnames(full.pheno) <- paste("Sim", 1:dim(full.pheno)[2], sep = "")
	
		#=================================================================	
		#add QTL
		#=================================================================
		if(num.qtl > 0){
			qtl.locale <- NULL
			num.genotypes <- apply(geno.obj$geno, 2, function(x) table(x))
			num.genotype <- matrix(NA, nrow = length(num.genotypes), ncol = max(unlist(lapply(num.genotypes, length))))
			for(i in 1:length(num.genotypes)){
				num.genotype[i,1:length(num.genotypes[[i]])] <- num.genotypes[[i]]
				}
			
			bad.loci <- which(num.genotype < 2, arr.ind = TRUE)
			good.loci <- setdiff(1:dim(num.genotype)[1], bad.loci[,1])
			
			chr <- geno.obj$chromosome[good.loci]
			u_chr <- unique(geno.obj$chromosome)
			if(num.qtl > length(u_chr)){replace = TRUE}else{replace = FALSE}
			qtl.chr <- sample(u_chr, num.qtl, replace = replace)
			if(num.qtl > num.pheno){replace = TRUE}else{replace = FALSE}
			qtl.pheno <- sample(1:num.pheno, num.qtl, replace = replace)
			rel.locale <- rep(NA, num.qtl)
		
			for(i in 1:length(qtl.chr)){
				chr.locale <- which(geno.obj$chromosome == qtl.chr[i])
				qtl.locale <- c(qtl.locale, chr.locale[sample(length(chr.locale), 1)])
				rel.locale[i] <- which(chr.locale == qtl.locale[i])
				genotypes <- unique(geno.obj$geno[,qtl.locale[i]])
				new.means <- segment.region(-10, 10, length(genotypes), "ends")
				for(g in 1:length(genotypes)){
					ind.g <- which(geno.obj$geno[,qtl.locale[i]] == genotypes[g])
					mean.pheno <- mean(full.pheno[ind.g,qtl.pheno[g]])
					sd.pheno <- sd(full.pheno[ind.g,qtl.pheno[g]])
					#change the means for the individuals with this genotype
					full.pheno[ind.g,qtl.pheno[g]] <- full.pheno[ind.g,qtl.pheno[g]] + rnorm(length(ind.g), new.means[g], sd.pheno)
					}
				#renormalize the phenotype	
				full.pheno[,qtl.pheno[g]] <- rz.transform(full.pheno[,qtl.pheno[g]])
				# quartz()
				# boxplot(full.pheno[,qtl.pheno[g]]~as.factor(geno.obj$geno[,qtl.locale]), main = paste("Pheno", qtl.pheno[i], " : Chr", qtl.chr[i]))
				table.order <- order(qtl.pheno)
				qtl.table <- cbind(colnames(full.pheno)[qtl.pheno[table.order]], qtl.chr[table.order], rel.locale[table.order])
				colnames(qtl.table) <- c("Pheno", "Chr", "Marker.Position")
				}
			}else{
				qtl.table <- "none"
				}

	get.means <- function(phen, geno){
		means <- unlist(lapply(levels(as.factor(geno)), function(x) mean(phen[which(geno == x)])))
		error <- unlist(lapply(levels(as.factor(geno)), function(x) sd(phen[which(geno == x)])))
		return(cbind(means,error))
		}

	get.vals <- function(phen, geno){
		geno.vals <- lapply(levels(as.factor(geno)), function(x) phen[which(geno == x)])
		return(geno.vals)
		}
		
	cols <- c("#7fc97f", "#beaed4", "#fdc086")
	# pdf("geno.effects.pdf", width = 30, height = 9)

	# for(ph in 1:num.pheno){
		# test <- apply(geno.obj$geno, 2, function(x) get.means(full.pheno[,ph], x))
	
		# ymax <- max(test[1:3,]+test[4:6,], na.rm = TRUE)
		# ymin <- min(test[1:3,]-test[4:6,], na.rm = TRUE)
		# a <- barplot(test[1:3,], beside = TRUE, border = FALSE, space = c(0, 2), col = cols, names = rep("", dim(geno.obj$geno)[2]), ylim = c(ymin, ymax), main = paste("Pheno", ph))
		# for(i in 1:3){
			# segments(x0 = a[i,], y0 = test[i,]-test[(i+3),], y1 = test[i,]+test[(i+3),], col = cols[i])
			# }
		# text(x = a[2,qtl.locale], y = rep(0, num.qtl), labels = "*", col = "red", cex = 2)
		# }
	# dev.off()
	
		# pdf("geno.effects.pdf", width = 30, height = 9)
		# for(ph in 1:num.pheno){
			# test <- apply(geno.obj$geno, 2, function(x) get.vals(full.pheno[,ph], x))		
			# plot.new()
			# plot.window(xlim = c(1,dim(geno.obj$geno)[2]*6), ylim = c(min(full.pheno[,ph], na.rm = TRUE), max(full.pheno[,ph], na.rm = TRUE)))
			# start.val <- 1
			# for(g in 1:length(test)){
				# at.vals <- segment.region(start.val+0.5, start.val+4-0.5, 3)
				# boxplot(test[[g]], at = at.vals, add = TRUE, col = cols, name = rep("", 3), axes = FALSE)
				# start.val = start.val + 6
				# if(length(which(qtl.locale == g)) > 0){
					# text(x = median(at.vals), y = min(full.pheno[,ph]), labels = "*", col = "red", cex = 2)
					# }
				# }
			# axis(2)
			# }
		# dev.off()
		
	rownames(full.pheno) <- 1:dim(full.pheno)[1]
	colnames(full.pheno) <- paste("P", 1:dim(full.pheno)[2], sep = "")
	
	
	return(list("pheno" = full.pheno, "QTL" = qtl.table))
	
	}