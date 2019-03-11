#This function compares the interaction terms from the pair wise regression
#with the direct influences


compare.interactions <- function(data.obj, p.or.q = 0.05){

	
	pairscan.results <- data.obj$pairscan.results
	dir.inf <- data.obj$var.to.var.influences
	
	pheno <- names(pairscan.results)
	marker.pairs <- pairscan.results[[1]][[1]][,1:2]
	
	interaction.stats <- list()
	for(ph in 1:length(pheno)){
		effect.table <- pairscan.results[[ph]][[1]]
		se.table <- pairscan.results[[ph]][[2]]
		effect <- as.numeric(effect.table[,dim(effect.table)[2]])
		se <- as.numeric(se.table[,dim(se.table)[[2]]])
		stat <- abs(effect/se)
		stat.table <- cbind(effect, se, stat)
		colnames(stat.table) <- c("coef", "se", "t")
		interaction.stats[[ph]] <- stat.table
		}
	names(interaction.stats) <- pheno
	
	#now gather the m12 and m21 stats in the same order as the interaction stats
	m12.stat.table <- NULL
	m21.stat.table <- NULL
	for(m in 1:length(marker.pairs[,1])){
		pair.locale <- intersect(which(dir.inf[,1] == marker.pairs[m,1]), which(dir.inf[,2] == marker.pairs[m,2]))
		m12.stats <- c(as.numeric(dir.inf[pair.locale,"m12"]), as.numeric(dir.inf[pair.locale,"m12.std.dev"]), as.numeric(dir.inf[pair.locale,"m12"])/as.numeric(dir.inf[pair.locale,"m12.std.dev"]))
		m12.stat.table <- rbind(m12.stat.table, m12.stats)
		
		m21.stats <- c(as.numeric(dir.inf[pair.locale,"m21"]), as.numeric(dir.inf[pair.locale,"m21.std.dev"]), as.numeric(dir.inf[pair.locale,"m21"])/as.numeric(dir.inf[pair.locale,"m21.std.dev"]))
		m21.stat.table <- rbind(m21.stat.table, m21.stats)
		
	}

	m.stats <- list(m12.stat.table, m21.stat.table)
	names(m.stats) <- c("m12", "m21")		

	max.coef <- apply(matrix(cbind(abs(m.stats[[1]][,1]), abs(m.stats[[2]][,1])), ncol = 2, byrow = FALSE), 1, max)
	max.stat <- apply(matrix(cbind(abs(m.stats[[1]][,3]), abs(m.stats[[2]][,3])), ncol = 2, byrow = FALSE), 1, max)

	#find the lines of significance for both sets.
	perm.pairscan.scan <- data.obj$pairscan.perm
	var.inf.sig <- data.obj$var.to.var.p.val
	
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.inf.sig))))

	
	dir.inf.sig.locale <- which(as.numeric(var.inf.sig[,var.sig.col]) <= p.or.q)
	if(length(dir.inf.sig.locale) > 0){
		m.sig.level <- min(as.numeric(var.inf.sig[dir.inf.sig.locale,"|Effect|/SE"]))
		}else{
		m.sig.level <- 10
		}
	
	get.pval <- function(stat, dist){
		pval <- length(which(abs(dist) >= abs(stat)))/length(dist)
		return(pval)
		}
	
	pairscan.scan.pvals <- list()
	for(ph in 1:length(pheno)){
		pairscan.scan.stat <- as.numeric(perm.pairscan.scan[[ph]][[1]][,dim(perm.pairscan.scan[[ph]][[1]])[2]])/as.numeric(perm.pairscan.scan[[ph]][[2]][,dim(perm.pairscan.scan[[ph]][[2]])[2]])
		pvals <- apply(matrix(interaction.stats[[ph]][,3], ncol = 1), 1, function(x) get.pval(x, pairscan.scan.stat))
		adj.pval <- p.adjust(pvals, method = "holm")
		sig.locale <- which(adj.pval <= p.or.q)
		if(length(sig.locale) > 0){
			pairscan.sig.level <- min(pairscan.scan.stat[adj.pval <= p.or.q])
			}else{
			pairscan.sig.level <- 10
			}
		pairscan.scan.pvals[[ph]] <- pairscan.sig.level
		}
	names(pairscan.scan.pvals) <- pheno



	#go through each phenotype
	layout.mat <- get.layout.mat((length(pheno)*3), type = "l")
	dev.new(width = dim(layout.mat)[2]*3, height = dim(layout.mat)[1]*3)
	# pdf(pdf.label, width = dim(layout.mat)[2]*3, height = dim(layout.mat)[1]*3)
	layout(layout.mat)
	for(ph in 1:length(pheno)){
		plot(interaction.stats[[ph]][,3], abs(m.stats[[1]][,3]), xlab = paste(pheno[ph], "stat"), ylab = "m12/se")
		abline(0,1); abline(h = m.sig.level, lty = 2); abline(v = pairscan.sig.level, lty = 2)
		if(ph == 1){
			par(xpd = TRUE)
			legend(0, max(abs(m.stats[[1]][,3]))*1.2, legend = "significance cutoff", lty = 2)
			par(xpd = FALSE)
			}
		
		plot(interaction.stats[[ph]][,3], abs(m.stats[[2]][,3]), xlab = paste(pheno[ph], "stat"), ylab = "m21/se")
		abline(0,1); abline(h = m.sig.level, lty = 2); abline(v = pairscan.sig.level, lty = 2)
		
		plot(interaction.stats[[ph]][,3], abs(max.stat), xlab = paste(pheno[ph], "stat"), ylab = "max(m12/se, m21/se)")
		abline(0,1); abline(h = m.sig.level, lty = 2); abline(v = pairscan.sig.level, lty = 2)
		}

}