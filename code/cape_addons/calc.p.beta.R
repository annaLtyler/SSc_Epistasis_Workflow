#This function calculates p values for each beta coefficient
#in the pairscan results. This is for comparing cape results
#to standard epistasis results. The function generates a 
#var.to.var.p.val table for each phenotype. If get.min.pval
#is set to TRUE, the function will merge the results for all
#phenotypes by finding the minimum interaction p value across
#all phenotoypes

calc.p.beta <- function(data.obj, pairscan.obj, pval.correction = c("holm", "fdr", "lfdr"), p.type = c("model", "empirical"), use.min.pval = FALSE, verbose = TRUE, run.parallel = TRUE, n.cores = 2){
	
	if(run.parallel){
		require(doParallel)
		require(parallel)
		require(foreach)
		}
	
	
	if(missing(pairscan.obj)){
		pairscan.results <- data.obj$pairscan.results
		pairscan.perm <- data.obj$pairscan.perm
		}else{
		pairscan.results <- pairscan.obj$pairscan.results
		pairscan.perm <- pairscan.obj$pairscan.perm
		}
	
	p.type.check <- grep("emp", p.type)
	if(length(p.type.check) > 0){p.type <- "empirical"}
		
	if(is.null(pairscan.perm[[1]]) && p.type == "empirical"){
		stop("Permutations must to calculate empirical p values")
		}
	
	get.p.val <- function(pair.val, null.dist){
		p.val <- length(which(abs(null.dist) >= abs(pair.val)))/length(null.dist)
		return(p.val)
		}
	
	
	pdf(paste("P_value_distributions_", p.type, ".pdf", sep = ""))
	all.var.to.var.p.val <- vector(mode = "list", length = length(pairscan.results))
	names(all.var.to.var.p.val) <- names(pairscan.results)
	for(i in 1:length(pairscan.results)){
		if(verbose){(cat("Calculating p values for ", names(pairscan.results)[i], "...\n", sep = ""))}
		last.pos <- dim(pairscan.results[[i]][[1]])[2]
		last.pos.perm <- dim(pairscan.perm[[i]][[1]])[2] 
		pair.vals <- matrix(as.numeric(pairscan.results[[i]][[1]][,last.pos])/as.numeric(pairscan.results[[i]][[2]][,last.pos]), ncol = 1)
		null.dist <- as.numeric(pairscan.perm[[i]][[1]][,last.pos.perm])/as.numeric(pairscan.perm[[i]][[2]][,last.pos.perm])
		if(p.type == "empirical"){
			if(run.parallel){
				# start.time <- Sys.time()
				cl <- makeCluster(n.cores)
				registerDoParallel(cl)
				pvals <- foreach(p = 1:length(pair.vals), .combine = "c") %dopar% {
				# pvals <- foreach(p = 1:1000, .combine = "c") %dopar% {
					get.p.val(pair.vals[p], null.dist)
					}
				stopCluster(cl)
				# stop.time <- Sys.time()
				# total.time <- stop.time - start.time
				# total.time
				}else{
				pvals <- apply(pair.vals, 1, function(x) get.p.val(x, null.dist))
				}
			}else{
			pvals <- pairscan.results[[i]][[4]][,3]
			}
			
		hist(pvals, main = names(pairscan.results)[i], breaks = 100)
		abline(v = 0.05, col = "red")

		var.to.var.p.val <- matrix(NA, nrow = nrow(pairscan.results[[i]][[1]]), ncol = 7)
		colnames(var.to.var.p.val) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical","p.adjusted")
		var.to.var.p.val[,c("Source", "Target")] <- pairscan.results[[i]][[1]][,1:2]
		var.to.var.p.val[,"Effect"] <- as.numeric(pairscan.results[[i]][[1]][,last.pos])
		var.to.var.p.val[,"SE"] <- as.numeric(pairscan.results[[i]][[2]][,last.pos])
		var.to.var.p.val[,"|Effect|/SE"] <- abs(as.numeric(pairscan.results[[i]][[1]][,last.pos])/as.numeric(pairscan.results[[i]][[2]][,last.pos]))
		var.to.var.p.val[,"P_empirical"] <- pvals

		if(length(pval.correction) > 1){
			pval.correction <- "holm"
			}

		if(pval.correction == "none"){
			var.to.var.p.val[,7] <- pvals
			colnames(var.to.var.p.val)[7] <- "p.adjusted"
			}
			
		if(pval.correction == "holm"){
			var.to.var.p.val[,7] <- p.adjust(pvals, method = "holm")
			colnames(var.to.var.p.val)[7] <- "p.adjusted"
			}
		
		if(pval.correction == "fdr" || pval.correction == "lfdr"){
			fdr.out <- fdrtool(pvals, statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "fndr")
			if(pval.correction == "lfdr"){
				lfdr <- fdr.out$lfdr
				var.to.var.p.val[,7] <- lfdr
				colnames(var.to.var.p.val)[7] <- "lfdr"
				}else{
				qval <- fdr.out$qval
				var.to.var.p.val[,7] <- qval
				colnames(var.to.var.p.val)[7] <- "qval"
				}
			}
		all.var.to.var.p.val[[i]] <- var.to.var.p.val
		}
	dev.off()
	#for each pair, find the phenotype in which its interaction
	#had the minimum p value consolidate these into a single table
	get.min.pval <- function(pair){
		pvals <- NULL
		for(i in 1:length(all.var.to.var.p.val)){
			pair.locale <- intersect(which(all.var.to.var.p.val[[i]][,1] == pair[1]), which(all.var.to.var.p.val[[i]][,2] == pair[2]))
			pvals <- c(pvals, all.var.to.var.p.val[[i]][pair.locale,7])
			}
		min.locale <- which(pvals == min(pvals, na.rm = TRUE))[1]
		return(all.var.to.var.p.val[[min.locale]][pair.locale,])
		}
	
	if(use.min.pval){
		if(verbose){cat("Finding minimum p value for each pair...\n")}
		all.pairs <- NULL
		for(i in 1:length(all.var.to.var.p.val)){
			all.pairs <- rbind(all.pairs, all.var.to.var.p.val[[i]][,1:2])
			}
		u_pairs <- unique(all.pairs)
		
		final.table <- t(apply(u_pairs, 1, get.min.pval))
	
		#replicate the table with the pairs in the opposite direction
		#because standard epistasis is symmetrical
		dup.table <- final.table
		dup.table[,1:2] <- final.table[,2:1]
		
		final.table <- rbind(final.table, dup.table)
		data.obj$var.to.var.p.val <- final.table
		}
	data.obj$all.var.to.var.p.val  <- all.var.to.var.p.val
	
	return(data.obj)
		
	
	
}