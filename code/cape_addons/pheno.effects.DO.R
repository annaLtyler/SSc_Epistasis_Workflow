#This function looks at the magnitude of phenotypic changes
#from main effects and from interactions (change from additivity)
# data.obj <- get.network(cross, p.or.q = 0.01, collapse.linked.markers = FALSE)

pheno.effects.DO <- function(data.obj, covar = NULL, scan.what = c("normalized.traits", "raw.traits"), geno.coding = c("Additive", "Dominant", "Recessive"), collapsed.net = FALSE, color.scheme = c("other", "DO/CC")){
	
	scan.what = scan.what[1]
	color.scheme = color.scheme[1]
	geno.coding = geno.coding[1]
	
	require(vioplot)
	if(is.null(data.obj$full.net)){
		stop("The full network must be present to run this script.")
		}
	#=========================================================
	#internal functions
	#=========================================================
	get.marker.name <- function(marker.name){
		if(length(grep("Chr", marker.name)) > 0){
			marker.locale <- which(names(data.obj$linkage.blocks.full) == marker.name)
			marker.label <- data.obj$linkage.blocks.full[[marker.locale]]
			return(marker.label)
			}else{
			if(is.null(allele)){
				full.name <- marker.name
				}else{
				full.name <- paste0(marker.name, "_", allele)
				}
			return(marker.name)
			}
		}


	get.motif.locale <- function(list1, list2){
		motif.list <- vector(mode = "list", length = length(list1))
		for(i in 1:length(list1)){
			motif.list[[i]] <- intersect(list1[[i]], list2[[i]])
			}
		return(motif.list)
		}

	get.markers <- function(block.name){
		block.locale <- which(names(data.obj$linkage.blocks.full) == block.name)
		return(data.obj$linkage.blocks.full[[block.locale]])	
		}
		
	
	
		get.motif.effects <- function(motif.names, motif.locale, collapsed.net){
			motif.effects <- vector(mode = "list", length = length(motif.locale))
			names(motif.effects) <- names(motif.names)
			for(i in 1:length(motif.locale)){
				if(length(motif.locale[[i]]) > 0){
					all.effects <- lapply(motif.locale[[i]], function(x) get.pheno.effects(data.obj, marker1.name = motif.names[[i]][x,1], marker2.name = motif.names[[i]][x,2], pheno.name = motif.names[[i]][x,3], scan.what = scan.what, collapsed.net = collapsed.net, geno.coding = geno.coding, covar = covar))
					effects.mat <- list2Matrix(all.effects)

					names.mat <- cbind(effects.mat[,1:2,drop=FALSE], rep(names(motif.names)[i], nrow(effects.mat)))
					effects.mat <- effects.mat[,3:5,drop=FALSE]
					colnames(effects.mat) <- c("m1", "m2", "interaction")
					rownames(effects.mat) <- apply(names.mat, 1, function(x) paste(x[1], x[2], x[3], sep = "|"))
					motif.effects[[i]] <- effects.mat
					}
				} #end looping through phenotypes
			return(motif.effects)
			}
	
	
		plot.motif.effects <- function(motif.effects, motif.name){
			plot.new()
			if(is.null(unlist(motif.effects))){
				maxy = 1; miny = 0
				}else{
				maxy <- max(unlist(motif.effects), na.rm = TRUE)	
				# miny <- min(unlist(motif.effects), na.rm = TRUE)
				miny <- 0
				}
			if(!is.finite(maxy)){maxy = 1; miny = 0}
			plot.window(xlim = c(0,1), ylim = c(0,1))
			text(0.5, 0.7, motif.name, cex = 2)
			text(x = c(0.2, 0.8), y = c(0.25, 0.25), labels = c("Add", "Int"), cex = 2)
			
			for(i in 1:length(motif.effects)){
					if(length(motif.effects[[i]]) > 0){
					add.int <- vector(mode = "list", length = 2)
					add.int[[1]] <- abs(motif.effects[[i]][,1] + motif.effects[[i]][,2])
					add.int[[2]] <- abs(motif.effects[[i]][,3])
					diff.results <- try(t.test(add.int[[1]], add.int[[2]]), silent = TRUE)
					if(class(diff.results) != "try-error"){
						avg.vals <- diff.results$estimate
						pval <- diff.results$p.value
						}else{
						avg.vals <- NA
						pval <- NA
						}
					par(bty = "n")
					plot.height <- maxy - miny
					plot.new()
					plot.window(xlim = c(0.7,2.3), ylim = c(miny, maxy))
					# plot(c(0:1), c(0:1), type = "n", xlim = c(0.7,2.3), ylim = c(0, plot.height*1.4), axes = FALSE, xlab = "", ylab = "")
					vioplot(add.int[[1]][which(!is.na(add.int[[1]]))], add.int[[2]][which(!is.na(add.int[[2]]))], col = "white", drawRect = FALSE, add = TRUE, at = c(1.2, 1.8), axes = FALSE, wex = 0.5)
					mtext(paste(names(motif.effects)[i], "\np =", signif(pval, 2)), line = 1)
					stripchart(list(add.int[[1]], add.int[[2]]), method = "jitter", add = TRUE, vertical = TRUE, pch = 16, at = c(1.2, 1.8))
					abline(h = 0)
					axis(2)
					par(xpd = TRUE)
					text(x = c(1.2, 1.8), y = rep(miny - (plot.height*0.15), 2), labels = c(signif(avg.vals[1], 2), signif(avg.vals[2],2)), cex = 2)
					}else{
					plot.new()
					plot.window(xlim = c(0,1), ylim = c(0,1))
					text(0.5, 0.5, "No motifs for\nthis phenotype")
					}
				} #end looping through phenotypes
				#show all main effects and interaction effects for all phenotypes
				all.pheno.add <- unlist(lapply(motif.effects, function(x) x[,1] + x[,2]))
				all.pheno.int <- unlist(lapply(motif.effects, function(x) x[,3]))
				if(length(all.pheno.add) > 0){
					vioplot(all.pheno.add[which(!is.na(all.pheno.add))], all.pheno.int, names = c("Main", "Int"), col = "white", drawRect = FALSE)
					mtext("All Phenotypes", cex = 1.5, line = 1)
					stripchart(list(all.pheno.add, all.pheno.int), method = "jitter", add = TRUE, vertical = TRUE, pch = 16)
					}else{
					plot.new()
					plot.window(xlim = c(0,1), ylim = c(0,1))
					text(0.5, 0.5, "No motifs for\nthis phenotype")
					}
			}
			
			
	
		#shuffle the edge weights in the network
		#keeping the topology the same
		shuffle.edges <- function(motif.net, sep.main){
			if(sep.main){
				gene.net.ind <- 1:(dim(motif.net)[1]-1)
				pheno.net.ind <- dim(motif.net)[1]
				gene.net <- motif.net[, gene.net.ind,drop=FALSE]
				pheno.net <- motif.net[, pheno.net.ind, drop=FALSE]
				int.locale <- which(gene.net != 0)
				int.weights <- gene.net[int.locale]
				main.locale <- which(pheno.net != 0)
				main.weights <- pheno.net[main.locale]
				gene.net[int.locale] <- sample(int.weights)
				pheno.net[main.locale] <- sample(main.weights)
				motif.net <- cbind(gene.net, pheno.net)
				}else{
				edge.locale <- which(motif.net != 0)
				edge.weights <- motif.net[edge.locale]
				motif.net[edge.locale] <- sample(edge.weights)
				}
			return(motif.net)
			}


		motif.mat <- function(motif.effects){
			row.total <- sum(unlist(lapply(motif.effects, nrow)))
			mat <- matrix(NA, nrow = row.total, ncol = 4)
			rownames(mat) <- rep("", row.total)
			block.start <- 1
			for(i in 1:length(motif.effects)){
				if(length(motif.effects[[i]]) > 0){
					row.ind <- block.start:(block.start+(nrow(motif.effects[[i]])-1))
					mat[row.ind,1:2] <- t(apply(motif.effects[[i]][,1:2,drop=FALSE], 1, function(x) sort(as.numeric(x), na.last = TRUE)))
					mat[row.ind,3] <- apply(motif.effects[[i]][,1:2,drop=FALSE], 1, function(x) sum(as.numeric(x), na.rm = TRUE))
					mat[row.ind,4] <- motif.effects[[i]][,3]
					
					rownames(mat)[row.ind] <- rownames(motif.effects[[i]])
					block.start <- block.start + nrow(motif.effects[[i]])
					}
				}
			return(mat)
			}


		plot.pheno.effect.lines <- function(motif.matrix, label){
			
			if(nrow(motif.matrix) == 0){
				plot.new()
				plot.window(xlim = c(0,1), ylim = c(0,1))
				text(0.5, 0.5, labels = paste("No", label, "Motifs"))
				return()
				}

			col1 <- get.color("blue")[2]; col2 <- get.color("brown")[2]; col3 <- "#de2d26"
			
			plot.height <- max(as.numeric(motif.matrix), na.rm = TRUE) - min(as.numeric(motif.matrix), na.rm = TRUE)
			colV <- rep(col1, nrow(motif.matrix))
			classV <- rep("less_than_additive", nrow(motif.matrix))
			
			int.act.diff <- abs(as.numeric(motif.matrix[,4])) - abs(as.numeric(motif.matrix[,3]))
			
			outside.add <- which(int.act.diff > 0)
			colV[outside.add] <- rep(col2, length(outside.add))
			classV[outside.add] <- rep("more_than_additive", length(outside.add))
			
			max.add <- max(abs(as.numeric(motif.matrix[,3])))
			extreme <- which(abs(as.numeric(motif.matrix[,4])) > max.add)
			colV[extreme] <- rep(col3, length(extreme))
			classV[extreme] <-  rep("extreme", length(extreme))
			 
			plot.new()
			plot.window(xlim = c(1,4), ylim = c(min(as.numeric(motif.matrix), na.rm = TRUE), max(as.numeric(motif.matrix), na.rm = TRUE)))
			for(i in 1:nrow(motif.matrix)){
				points(as.numeric(motif.matrix[i,]), type = "b", col = colV[i], pch = c(1,1,1,16))
				}
			# apply(motif.matrix, 1, function(x) points(x, type = "b"))
			axis(2)
			mtext(label, cex = 2)
			par(xpd = TRUE)
			text(x = c(1:4), y = rep(min(as.numeric(motif.matrix), na.rm = TRUE)-(plot.height*0.1), 4), labels = c("Main1", "Main2", "Additive", "Actual"), cex = 1.7)
			par(xpd = FALSE)
			
			abline(h = 0, lty = 2)
			
			num.motifs <- nrow(motif.matrix)
			
			class.count <- table(classV)
			cat(label, "\n")
			
			cat(paste(names(class.count), signif(class.count/sum(class.count), 2), "\t"), "\n")
			
			marker.info <- strsplit(rownames(motif.matrix), "\\|")
			source.marker <- sapply(marker.info, function(x) x[1])
			target.marker <- sapply(marker.info, function(x) x[2])
			target.pheno <- sapply(marker.info, function(x) x[3])
			matrix.with.classes <- cbind(source.marker, target.marker, target.pheno, motif.matrix, classV)
			colnames(matrix.with.classes) <- c("Source", "Target", "Pheno", "Main1", "Main2", "Additive", "Actual", "Class")
			return(matrix.with.classes)
			}
			
		report.stats <- function(labeled.motif.mat){
			layout.matrix <- matrix(c(1:4, 4, 5), nrow = 2, byrow = TRUE)
			layout(layout.matrix)
	
			#look at overall stats
			phenos.affected <- table(labeled.motif.mat[,"Pheno"])
			barplot(phenos.affected, main = "Overall Phenotypes Affected", cex.axis = 2)
	
			source.alleles <- sapply(strsplit(labeled.motif.mat[,"Source"], "_"), function(x) x[2])
			source.table <- table(source.alleles)
			barplot(source.table, col = allele.colors[match(names(source.table), allele.colors[,2]), 3], main = "Overall Source Alleles", names = allele.colors[match(names(source.table), allele.colors[,2]), 1], cex.names = 1.5, cex.axis = 2)
	
			target.alleles <- sapply(strsplit(labeled.motif.mat[,"Target"], "_"), function(x) x[2])
			target.table <- table(target.alleles)
			barplot(target.table, col = allele.colors[match(names(target.table), allele.colors[,2]), 3], main = "Overall Target Alleles", names = allele.colors[match(names(target.table), allele.colors[,2]), 1], cex.names = 1.5, cex.axis = 2)
	
			allele.combos <- apply(cbind(source.alleles, target.alleles), 1, function(x) paste(x[1], x[2], sep = "_"))
			combo.table <- table(allele.combos)
			barplot(combo.table, main = "Allele Combinations (Source_Target)")
					
			diff.locale <- apply(cbind(source.alleles, target.alleles), 1, function(x) x[1] != x[2])
			frac.diff <- length(which(diff.locale))/length(source.alleles)
			barplot(frac.diff, main = paste("Proportion Different Parents\n", signif(frac.diff, 2)), ylim = c(0, 1))
				
			#and stats for individual classes
			classes <- unique(labeled.motif.mat[,"Class"])
			classes <- c(classes, "all_more_than_additive")
			class.results <- vector(mode = "list", length = length(classes))
			names(class.results) <- classes
			for(i in 1:length(classes)){
				layout(layout.matrix)
				
				if(classes[i] == "all_more_than_additive"){
					class.locale <- c(which(labeled.motif.mat[,"Class"] == "more_than_additive"), which(labeled.motif.mat[,"Class"] == "extreme"))	
					}else{
					class.locale <- which(labeled.motif.mat[,"Class"] == classes[i])	
					}
					
				class.mat <- labeled.motif.mat[class.locale,,drop=FALSE]
				phenos.affected <- table(class.mat[,"Pheno"])
				barplot(phenos.affected, main = paste(classes[i], "\nPhenotypes Affected"), cex.axis = 2)
				
				source.alleles <- sapply(strsplit(class.mat[,"Source"], "_"), function(x) x[2])
				source.table <- table(source.alleles)
				barplot(source.table, col = allele.colors[match(names(source.table), allele.colors[,2]), 3], main = paste(classes[i], "\nSource Alleles"), names = allele.colors[match(names(source.table), allele.colors[,2]), 1], cex.axis = 2, cex.names = 1.5)
				
				target.alleles <- sapply(strsplit(class.mat[,"Target"], "_"), function(x) x[2])
				target.table <- table(target.alleles)
				barplot(target.table, col = allele.colors[match(names(target.table), allele.colors[,2]), 3], main = paste(classes[i], "\nTarget Alleles"), names = allele.colors[match(names(target.table), allele.colors[,2]), 1], cex.names = 1.5, cex.axis = 2)

				allele.combos <- apply(cbind(source.alleles, target.alleles), 1, function(x) paste(x[1], x[2], sep = "_"))
				combo.table <- table(allele.combos)
				barplot(combo.table, main = "Allele Combinations (Source_Target)")
				
				diff.locale <- apply(cbind(source.alleles, target.alleles), 1, function(x) x[1] != x[2])
				frac.diff <- length(which(diff.locale))/length(source.alleles)
				barplot(frac.diff, main = paste("Proportion Different Parents\n", signif(frac.diff, 2)), ylim = c(0, 1))
	
				}
			
			return(class.results)
			}
		#================================================================
		allele.colors <- get.allele.colors(color.scheme = color.scheme)
		pheno <- get.pheno(data.obj, scan.what = scan.what, covar = covar)
		
		covar.info <- get.covar(data.obj)
		if(!is.null(covar.info$covar.table)){
			not.na.locale <- which(!is.na(rowSums(covar.info$covar.table)))
			}else{
			not.na.locale <- 1:nrow(pheno)	
			}

		#find all the motifs so we can break up the main effects
		#and interaction effects by motif type	
		motifs <- find.motifs(data.obj, collapsed.net = collapsed.net)
		geno <- data.obj$geno.for.pairscan

		motif.dir <- motifs[[2]]
		
		enhancing.locale <- lapply(motif.dir, function(x) which(x[,1] == 1))		
		suppressing.locale <- lapply(motif.dir, function(x) which(x[,1] == -1))
		coherent.locale <- lapply(motif.dir, function(x) which(x[,2] == x[,3]))
		incoherent.locale <- lapply(motif.dir, function(x) which(x[,2] != x[,3]))
		
		en.coh <- get.motif.locale(enhancing.locale, coherent.locale)
		en.inc <- get.motif.locale(enhancing.locale, incoherent.locale)
		supp.coh <- get.motif.locale(suppressing.locale, coherent.locale)
		supp.inc <- get.motif.locale(suppressing.locale, incoherent.locale)
	

		en.coh.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = en.coh, collapsed.net)
		en.inc.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = en.inc, collapsed.net)
		supp.coh.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = supp.coh, collapsed.net)
		supp.inc.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = supp.inc, collapsed.net)

		
		en.coh.mat <- motif.mat(en.coh.effects)
		en.inc.mat <- motif.mat(motif.effects = en.inc.effects)
		supp.coh.mat <- motif.mat(supp.coh.effects)		
		supp.inc.mat <- motif.mat(supp.inc.effects)
		
				
		pdf(paste0("Motif.Effects.Lines.", geno.coding, ".pdf"), width = 10, height = 10)
		par(mfrow = c(2,2))
		en.coh.mat.c <- plot.pheno.effect.lines(motif.matrix = en.coh.mat, "Enhancing Coherent")
		en.inc.mat.c <- plot.pheno.effect.lines(motif.matrix = en.inc.mat, "Enhancing Incoherent")
		supp.coh.mat.c <- plot.pheno.effect.lines(supp.coh.mat, "Suppressing Coherent")
		supp.inc.mat.c <- plot.pheno.effect.lines(supp.inc.mat, "Suppressing Incoherent")
		dev.off()
		
		write.table(en.inc.mat.c[order(as.numeric(en.inc.mat.c[,4])),,drop=FALSE], file = paste0("Motifs.Enhancing.Incoherent.", geno.coding, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(en.coh.mat.c[order(as.numeric(en.coh.mat.c[,4])),,drop=FALSE], file = paste0("Motifs.Enhancing.Coherent.", geno.coding, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(supp.coh.mat.c[order(as.numeric(supp.coh.mat.c[,4])),,drop=FALSE], file = paste0("Motifs.Suppressing.Coherent.", geno.coding, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(supp.inc.mat.c[order(as.numeric(supp.inc.mat.c[,4])),,drop=FALSE], file = paste0("Motifs.Suppressing.Incoherent.", geno.coding, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	final.result <- list(en.coh.mat.c, en.inc.mat.c, supp.coh.mat.c, supp.inc.mat.c)
	names(final.result) <- c("en.coh", "en.inc", "supp.coh", "supp.inc")
	invisible(final.result)

		
}