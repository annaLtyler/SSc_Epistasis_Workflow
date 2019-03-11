#This function looks at epistasis in animals with extreme phenotypes
# setwd("~/Documents/Data/Little_Cross/little_cross_data_cape/Body_Comp_IGF_New_P")
# setwd("~/Documents/Data/Frankel/Combined_Cross/Results/Combined_Cross_New_P_no_thresh")
# data.obj <- readRDS("cross.RData"); geno.obj <- NULL
# rownames(data.obj$pheno) <- rownames(data.obj$geno) <- 1:nrow(data.obj$pheno)
# data.obj <- marker2covar(data.obj, markers = c("sex", "Final.IGF.1"))
# data.obj <- get.network(data.obj, p.or.q = 0.01, standardized = FALSE, collapse.linked.markers = FALSE, verbose = TRUE)

extreme.pheno <- function(data.obj, geno.obj = NULL, min.pheno.percentile = 90, max.pheno.percentile = 100, main.effect1 = 1, main.effect2 = 1, interaction.effect = 1, alt.val = 0.5, nperm = 1000){

	
	pheno <- data.obj$raw.pheno
	geno <- get.geno(data.obj, geno.obj)
	num.pheno <- dim(pheno)[2]
		
	motifs <- find.motifs(data.obj, collapsed.net = FALSE, include.covar = FALSE)	


	syn <- lapply(motifs[[2]], function(x) Reduce("intersect", list(which(x[,1] == main.effect1), which(x[,2] == main.effect2), which(x[,3] == interaction.effect))))
	syn2 <- lapply(motifs[[2]], function(x) Reduce("intersect", list(which(x[,1] == main.effect2), which(x[,2] == main.effect1), which(x[,3] == interaction.effect))))

	for(i in 1:length(syn)){
		syn[[i]] <- unique(c(syn[[i]], syn2[[i]]))
		}


	min.val <- apply(pheno, 2, function(x) get.percentile(x, min.pheno.percentile))
	max.val <- apply(pheno, 2, function(x) get.percentile(x, max.pheno.percentile))
	ind <- apply(pheno, 2, function(x) intersect(which(x >= min.val), which(x <= max.val)))


	
	#====================================================================================
	# internal functions
	#====================================================================================
	get.genotypes <- function(ind.idx, marker1.idx, marker2.idx){
		all.ind.genos <- vector(mode = "list", length = length(ind.idx))
		for(i in 1:length(ind.idx)){
			ind.genos <- matrix(NA, nrow = length(marker1.idx), ncol = 2)
			for(j in 1:length(marker1.idx)){
				ind.genos[j,] <- c(geno[ind.idx[i], marker1.idx[j]], geno[ind.idx[i], marker2.idx[j]])
				}
			all.ind.genos[[i]] <- ind.genos
			}
		return(all.ind.genos)	
		}
		
	#alt.val denotes probability of alternate allele at that locus
	#this function finds the proportion of locus pairs for which each
	#individual has the alternate allele at both loci
	tally.alt.allele <- function(genotypes, alt.val){
		num.loci <- unlist(lapply(genotypes, nrow))
		alt.locus <- lapply(genotypes, function(x) x > alt.val)
		both.alt <- unlist(lapply(alt.locus, function(x) length(which(rowSums(x) == 2))))
		prop.both.alt <- both.alt/num.loci
		return(prop.both.alt)
		}

	process.motif.list <- function(motif.list, ind.which){
		if(nrow(motif.list) == 0){
			return(NULL)
			}
		source.idx <- get.marker.idx(data.obj, unlist(data.obj$linkage.blocks.full[motif.list[,1]]))
		target.idx <- get.marker.idx(data.obj, unlist(data.obj$linkage.blocks.full[motif.list[,2]]))
			
		#Find the proportion of synergistic interactions for which
		#these individuals have the alternate allele at both loci
		ind.genos <- get.genotypes(ind.which, source.idx, target.idx)
		prop.both.alt <- tally.alt.allele(ind.genos, alt.val)
		return(prop.both.alt)
		}


	permute.motif.list <- function(motif.list, ind.which){
		if(nrow(motif.list) == 0){
			return(NULL)
			}
		source.idx <- get.marker.idx(data.obj, unlist(data.obj$linkage.blocks.full[motif.list[,1]]))
		target.idx <- get.marker.idx(data.obj, unlist(data.obj$linkage.blocks.full[motif.list[,2]]))
			
		prop.perm.by.ind <- matrix(nrow = length(ind.which), ncol = nperm)
		cat("\n")
		for(n in 1:nperm){
			report.progress(n, nperm)
			rnd.ind <- sample(1:nrow(pheno), length(ind.which))
			rnd.genos <- get.genotypes(rnd.ind, source.idx, target.idx)
			prop.perm.by.ind[,n] <- tally.alt.allele(rnd.genos, alt.val)
			}
		return(prop.perm.by.ind)
		}
		
	plot.results <- function(motif.results, perm.results){
		motif.means <- lapply(motif.results, function(x) if(length(x) > 0){mean(x)})
		layout.mat <- get.layout.mat(length(motif.results))
		layout(layout.mat)
		for(i in 1:length(motif.results)){
			if(length(perm.results[[i]]) > 0){
				perm.means <- apply(perm.results[[i]], 2, mean)
				hist(perm.means, xlim = c(0, 0.5), main = colnames(pheno)[i], breaks = 100, xlab = "Proportion Pairs with Alternate")
				abline(v = motif.means[[i]], col = "red")	
				# mtext(label, outer = TRUE, line = -1.5)
				}
			}
		}
	#====================================================================================


	#high phenotype, synergistic positive interactions
	prop.motif.alt <- vector(mode = "list", length = ncol(pheno))
	prop.perm <- vector(mode = "list", length = ncol(pheno))
	names(prop.motif.alt) <- names(prop.perm) <- colnames(pheno)
	for(i in 1:ncol(pheno)){
		cat("\n", colnames(pheno)[i], "\n")
		prop.motif.alt[[i]] <- process.motif.list(motif.list = motifs[[1]][[i]][syn[[i]],,drop=FALSE], ind.which = ind[[i]])
		prop.perm[[i]] <- permute.motif.list(motifs[[1]][[i]][syn[[i]],,drop=FALSE], ind[[i]])
		}
	plot.results(prop.motif.alt, prop.perm)

		
	}