#This function looks at the magnitude of phenotypic changes
#from main effects and from interactions (change from additivity)
# data.obj <- get.network(cross, p.or.q = 0.01, collapse.linked.markers = FALSE)
#it is not fully migrated over from the DO version. the full network works
#as far as plot.lines

pheno.effects <- function(data.obj, covar = NULL, collapsed.net = TRUE, include.covar = FALSE, significant.main.effects.only = TRUE){
	
	require(vioplot)
	
	if(!significant.main.effects.only){
		new.data.obj <- get.network(data.obj, 1, collapse.linked.markers = collapsed.net)
		if(collapsed.net){
			filled.net <- new.data.obj$collapsed.net
			full.main.effects <- filled.net[,(nrow(filled.net)+1):ncol(filled.net)]
			new.net <- data.obj$collapsed.net
			new.net[,(nrow(new.net)+1):ncol(new.net)] <- full.main.effects
			data.obj$collapsed.net <- new.net
			}else{
			filled.net <- new.data.obj$full.net
			full.main.effects <- filled.net[,(nrow(filled.net)+1):ncol(filled.net)]
			new.net <- data.obj$full.net
			new.net[,(nrow(new.net)+1):ncol(new.net)] <- full.main.effects
			data.obj$full.net <- new.net
			}	
		}
	
	#=========================================================
	#internal functions
	#=========================================================
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
		
	
	get.effects <- function(m1, m2, pheno.name, min.hom = 0.5, collapsed.net){
		if(collapsed.net){
			m1.markers <- data.obj$linkage.blocks.collapsed[[which(names(data.obj$linkage.blocks.collapsed) == m1)]]
			m2.markers <- data.obj$linkage.blocks.collapsed[[which(names(data.obj$linkage.blocks.collapsed) == m2)]]
			m1.ind <- match(m1.markers, colnames(geno))
			m2.ind <- match(m2.markers, colnames(geno))
			all.int <- data.obj$full.net[m1.ind, m2.ind,drop=FALSE]
			max.int <- which(abs(all.int) == max(abs(all.int)), arr.ind = TRUE)
			m1.geno <- geno[,m1.ind[max.int[1]]]
			m2.geno <- geno[,m2.ind[max.int[2]]]	
			}else{
			m1.ind <- which(names(data.obj$linkage.blocks.full) == m1)
			m2.ind <- which(names(data.obj$linkage.blocks.full) == m2)
			m1.geno <- geno[,m1.ind]
			m2.geno <- geno[,m2.ind]
			}

		if(ncol(data.obj$p.covar.table) > 1){
			not.na.list <- apply(data.obj$p.covar.table, 2, function(x) which(!is.na(x)))
			if(class(not.na.list) == "matrix"){
			not.na.locale <- not.na.list[,1]				
			}else{
			not.na.locale <- Reduce("intersect", apply(data.obj$p.covar.table, 2, function(x) which(!is.na(x))))
			}
		}else{
		not.na.locale <- which(!is.na(data.obj$p.covar.table[,1]))	
		}
		

		pheno.locale <- which(colnames(pheno) == pheno.name)
		# pheno.vals <- pheno[,pheno.locale]
	
		add1.model <- lm(pheno[,pheno.locale]~m1.geno[not.na.locale])
		add2.model <- lm(pheno[,pheno.locale]~m2.geno[not.na.locale])
		m1.coef <- coefficients(add1.model)[2]
		m2.coef <- coefficients(add2.model)[2]
		add.prediction <- m1.coef + m2.coef
		
		int.model <- lm(pheno[,pheno.locale]~m1.geno[not.na.locale]*m2.geno[not.na.locale])
		actual.int <- sum(coefficients(int.model)[2:4])
		
		# m1.alt.hom <- pheno.vals[which(m1.geno >= min.hom)]
		# m1.ref.hom <- pheno.vals[which(m1.geno < min.hom)]
		# m2.alt.hom <- pheno.vals[which(m2.geno >= min.hom)]
		# m2.ref.hom <- pheno.vals[which(m2.geno < min.hom)]

		# # boxplot(list(c(m1.ref.hom, m2.ref.hom), m1.alt.hom, m2.alt.hom))

		# #find the difference in phenotype between animals with the alt allele
		# #and animals with the ref allele (all other alleles)
		# #This is the main effect of individual allele
		# m1.effect <- mean(m1.alt.hom, na.rm = TRUE) - mean(m1.ref.hom, na.rm = TRUE)
		# m2.effect <- mean(m2.alt.hom, na.rm = TRUE) - mean(m2.ref.hom, na.rm = TRUE)
		# add.prediction <- m1.effect + m2.effect

		#find the effect of the interaction
		# m1.m2.ref.hom <- pheno.vals[intersect(which(m1.geno < min.hom), which(m2.geno < min.hom))]
		# m1.m2.alt.hom <- pheno.vals[intersect(which(m1.geno >= min.hom), which(m2.geno >= min.hom))]
		# actual.int <- mean(m1.m2.alt.hom, na.rm = TRUE) - mean(m1.m2.ref.hom, na.rm = TRUE)
						
		return(c(m1.coef, m2.coef, actual.int))
		}


	
		get.motif.effects <- function(motif.names, motif.locale, min.hom, collapsed.net){
			motif.effects <- vector(mode = "list", length = length(motif.locale))
			names(motif.effects) <- names(motif.names)
			for(i in 1:length(motif.locale)){
				if(length(motif.locale[[i]]) > 0){
					all.effects <- lapply(motif.locale[[i]], function(x) get.effects(m1 = motif.names[[i]][x,1], m2 = motif.names[[i]][x,2], pheno.name = motif.names[[i]][x,3], min.hom, collapsed.net))
					motif.effects[[i]] <- list2Matrix(all.effects)
					colnames(motif.effects[[i]]) <- c("m1", "m2", "interaction")
					rownames(motif.effects[[i]]) <- apply(motif.names[[i]][motif.locale[[i]],,drop=FALSE], 1, function(x) paste(x[1], x[2], x[3], sep = "|"))
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
					mat[row.ind,1:2] <- t(apply(motif.effects[[i]][,1:2,drop=FALSE], 1, function(x) sort(x, na.last = TRUE)))
					mat[row.ind,3] <- apply(motif.effects[[i]][,1:2,drop=FALSE], 1, function(x) sum(x, na.rm = TRUE))
					mat[row.ind,4] <- motif.effects[[i]][,3]
					
					rownames(mat)[row.ind] <- rownames(motif.effects[[i]])
					block.start <- block.start + nrow(motif.effects[[i]])
					}
				}
			return(mat)
			}


		plot.lines <- function(motif.matrix, label){

			col1 <- get.color("blue")[2]; col2 <- get.color("brown")[2]; col3 <- "#de2d26"
			
			plot.height <- max(motif.matrix, na.rm = TRUE) - min(motif.matrix, na.rm = TRUE)
			colV <- rep(col1, nrow(motif.matrix))
			classV <- rep("less_than_additive", nrow(motif.matrix))
			
			outside.add <- which(apply(motif.matrix, 1, function(x) abs(x[3]) < abs(x[4])))
			colV[outside.add] <- rep(col2, length(outside.add))
			classV[outside.add] <- rep("more_than_additive", length(outside.add))
			
			extreme <- c(which(motif.matrix[,4] < min(motif.matrix[,3])), which(motif.matrix[,4] > max(motif.matrix[,3])))
			colV[extreme] <- rep(col3, length(extreme))
			classV[extreme] <-  rep("extreme", length(extreme))
			 
			plot.new()
			plot.window(xlim = c(1,4), ylim = c(min(motif.matrix, na.rm = TRUE), max(motif.matrix, na.rm = TRUE)))
			for(i in 1:nrow(motif.matrix)){
				points(motif.matrix[i,], type = "b", col = colV[i], pch = c(1,1,1,16))
				}
			# apply(motif.matrix, 1, function(x) points(x, type = "b"))
			axis(2)
			mtext(label, cex = 2)
			par(xpd = TRUE)
			text(x = c(1:4), y = rep(min(motif.matrix, na.rm = TRUE)-(plot.height*0.1), 4), labels = c("Main1", "Main2", "Additive", "Actual"), cex = 1.7)
			par(xpd = FALSE)
			
			num.motifs <- nrow(motif.matrix)
			
			class.count <- table(classV)
			print(class.count/sum(class.count))
			
			marker.info <- strsplit(rownames(motif.matrix), "\\|")
			source.marker <- sapply(marker.info, function(x) x[1])
			target.marker <- sapply(marker.info, function(x) x[2])
			target.pheno <- sapply(marker.info, function(x) x[3])
			matrix.with.classes <- cbind(source.marker, target.marker, target.pheno, motif.matrix, classV)
			colnames(matrix.with.classes) <- c("Source", "Target", "Pheno", "Main1", "Main2", "Additive", "Actual", "Class")
			return(matrix.with.classes)
			}
			
	
		#================================================================


		pheno <- get.pheno(data.obj, scan.what = "raw.traits", covar)

		#find all the motifs so we can break up the main effects
		#and interaction effects by motif type	
		motifs <- find.motifs(data.obj, collapsed.net = collapsed.net, include.covar = include.covar)
		
		if(sum(unlist(lapply(motifs[[1]], nrow))) == 0){stop("No motifs found.")}
		
		covar.table <- data.obj$p.covar.table
		colnames(covar.table) <- data.obj$p.covar
		geno <- cbind(data.obj$geno.for.pairscan, covar.table)
		
		motif.dir <- motifs[[2]]
		
		enhancing.locale <- lapply(motif.dir, function(x) which(x[,1] == 1))		
		suppressing.locale <- lapply(motif.dir, function(x) which(x[,1] == -1))
		coherent.locale <- lapply(motif.dir, function(x) which(x[,2] == x[,3]))
		incoherent.locale <- lapply(motif.dir, function(x) which(x[,2] != x[,3]))
		
		en.coh <- get.motif.locale(enhancing.locale, coherent.locale)
		en.inc <- get.motif.locale(enhancing.locale, incoherent.locale)
		supp.coh <- get.motif.locale(suppressing.locale, coherent.locale)
		supp.inc <- get.motif.locale(suppressing.locale, incoherent.locale)
	
		en.coh.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = en.coh, min.hom, collapsed.net)
		en.inc.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = en.inc, min.hom, collapsed.net)
		supp.coh.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = supp.coh, min.hom, collapsed.net)
		supp.inc.effects <- get.motif.effects(motif.names = motifs[[1]], motif.locale = supp.inc, min.hom, collapsed.net)
		
		en.coh.mat <- motif.mat(en.coh.effects)
		en.inc.mat <- motif.mat(motif.effects = en.inc.effects)
		supp.coh.mat <- motif.mat(supp.coh.effects)		
		supp.inc.mat <- motif.mat(supp.inc.effects)
		
				
		pdf("pheno.effects.pdf", width = 9, height = 6)
		par(mfrow = c(2,2))
		en.coh.mat.c <- plot.lines(motif.matrix = en.coh.mat, "Enhancing Coherent")
		en.inc.mat.c <- plot.lines(motif.matrix = en.inc.mat, "Enhancing Incoherent")
		supp.coh.mat.c <- plot.lines(supp.coh.mat, "Suppressing Coherent")
		supp.inc.mat.c <- plot.lines(supp.inc.mat, "Suppressing Incoherent")
		dev.off()
		
		
		write.table(en.coh.mat.c[order(as.numeric(en.coh.mat.c[,4])),,drop=FALSE], file = "Motifs.Enhancing.Coherent.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(en.inc.mat.c[order(as.numeric(en.inc.mat.c[,4])),,drop=FALSE], file = "Motifs.Enhancing.Incoherent.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(supp.coh.mat.c[order(as.numeric(supp.coh.mat.c[,4])),,drop=FALSE], file = "Motifs.Suppressing.Coherent.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(supp.inc.mat.c[order(as.numeric(supp.inc.mat.c[,4])),,drop=FALSE], file = "Motifs.Suppressing.Incoherent.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		
	
		
}