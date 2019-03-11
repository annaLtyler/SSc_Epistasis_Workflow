#This function shows uncollapsed main effects the way plotNetwork does, 
#but shows interactions in a different way. It's best if this just plots
#one chromosome at a time. It draws a colored line underneath the 
#interacting markers to indicate which other chromosome they interact 
#with. This is mostly for judging how well main effects and interaction
#effects line up in the blocks

plotNetwork_complex <- function(data.obj, chr = 1, interaction.line.width = 1){

	# require(igraph)

	if(length(chr) > 1){
		stop("This script plots only one chromosome at a time.")
		}

	collapsed.net = FALSE
	show.effect.size = FALSE
	trait = NULL
			
	all.chr <- data.obj$chromosome
	all.pos <- data.obj$marker.location

	adj.mat <- data.obj$full.net	

	if(is.null(adj.mat)){
		stop("get.network() must be run before plotting.")
		}
		
	blocks <- data.obj$linkage.blocks.full	
			
	all.markers <- as.vector(unlist(blocks))
		
		
		#assign a chromosome and relative position to each block
		get.chr.pos <- function(block){
			marker.locale <- which(colnames(data.obj$geno) %in% block)
			chr <- unique(all.chr[marker.locale])
			if(length(chr) > 1){
				chr.char <- paste(chr, collapse = ", ")
				stop(paste("There is linkage between markers on chromosomes ", chr.char,". Please try a high r2.thresh.", sep = ""))
				}
			pos <- mean(all.pos[marker.locale])
			max.pos <- max(all.pos[all.chr == chr])
			return(c(chr, pos/max.pos))
			}
		
		chr.pos <- t(sapply(blocks, get.chr.pos))
		colnames(chr.pos) <- c("chromosome", "position")

		used.chr <- as.numeric(chr.pos[,"chromosome"])

		if(is.null(chr)){
			chr <- sort(unique(used.chr))
			}else{
			chr <- sort(chr)	
			}
			
		#calculate beginning and end x coordinates for each chromosome
		num.chr <- length(chr)
		chr.x <- matrix(NA, ncol = 2, nrow = num.chr)
		rownames(chr.x) <- chr
		chr.x[,1] <- (0:(num.chr-1))+(num.chr*0.005)
		chr.x[,2] <- (1:num.chr)-(num.chr*0.005)


		#get the average effect size for the var to 
		#pheno effects for each block
		var.to.pheno <- data.obj$max.var.to.pheno.influence

		get.block.effect <- function(block){
			effects <- rep(NA, length(var.to.pheno)); names(effects) <- names(var.to.pheno)
			marker.locale <- lapply(var.to.pheno, function(x) match(block, x[,1]))
			for(i in 1:length(marker.locale)){
				effects[i] <- mean(var.to.pheno[[i]][marker.locale[[i]], "|t.stat|"])
				}
			return(effects)	
			}
			
		block.effects <- t(sapply(blocks, get.block.effect))
		rel.block.effects <- block.effects/max(block.effects)
		
		#and each phenotype
		if(is.null(trait)){
			pheno <- names(data.obj$max.var.to.pheno.influence)
			}else{
			pheno <- trait
			trait.locale <- which(trait %in% names(data.obj$max.var.to.pheno.influence))
			if(length(trait.locale) < length(trait)){
				not.found <- which(!trait %in% names(data.obj$max.var.to.pheno.influence))
				message("I couldn't find the following traits:")
				cat(trait[not.found], sep = "\n")
				return()
				}
			}
			
		
		#if we need to filter chr.pos and adj.mat to include
		#only the chromomsomes and phenotypes we are including
		chr.pos <- chr.pos[which(chr.pos[,"chromosome"] %in% chr),,drop = FALSE]
		
		if(dim(chr.pos)[1] > 1){
		all.block.pheno <- c(rownames(chr.pos), pheno)
		
		#keep the full adjacency matrix, since we'll need it to find all the interactions
		# adj.mat <- adj.mat[,colnames(adj.mat) %in% all.block.pheno, drop = FALSE]
		# adj.mat <- adj.mat[rownames(adj.mat) %in% rownames(chr.pos),, drop = FALSE]
		
		#get the absolute position of each marker
		get.abs.mrk <- function(chr.pos.row){
			chr.locale <- which(rownames(chr.x) == chr.pos.row[1])
			mrk.pos <- ((chr.x[chr.locale,2] - chr.x[chr.locale,1])*as.numeric(chr.pos.row[2])) + chr.x[chr.locale,1]
			return(as.numeric(mrk.pos))
			}
		
		ph.x <- apply(chr.pos, 1, get.abs.mrk)
		
		#if there are covariates, expand the covariate "chromosome"
		#to give each covariate its own segment
		covar.exists <- which(chr.pos[,1] == 0)
		covar.locale <- which(rownames(chr.x) == 0)
		covar.names <- names(covar.exists)
		if(length(covar.locale) > 0){
			covar.min <- chr.x[covar.locale,1]
			covar.max <- chr.x[covar.locale,2]
			new.cov.x <- matrix(NA, ncol = 2, nrow = length(which(used.chr == 0)))
			region.centers <- segment.region(covar.min, covar.max, length(which(used.chr == 0)), "center") 
			new.cov.x[,1] <- region.centers - ((covar.max-covar.min)/100)
			new.cov.x[,2] <- region.centers + ((covar.max-covar.min)/100)
			rownames(new.cov.x) <- rep(0, dim(new.cov.x)[1])
			
			#replace the original covar chromosome location with the new chromosome locations
			chr.x <- chr.x[which(rownames(chr.x) != 0),,drop=FALSE] #remove the covariate chromosomes from chr.x
			chr.x <- rbind(chr.x, new.cov.x) #put the covariate chromosome locations at the end of chr.x
			
			#also adjust the positions of the covariates in ph.x
			covar.blocks <- rownames(chr.pos)[which(chr.pos[,1] == 0)]
			if(length(covar.blocks) > 0){
				sig.covar <- NULL
				for(i in 1:length(covar.blocks)){
					sig.covar <- c(sig.covar, blocks[[covar.blocks[i]]])
					}
				# covar.pos <- which(colnames(data.obj$geno)[which(used.chr == 0)] %in% sig.covar)
				cov.centers <- segment.region(covar.min, covar.max, length(which(used.chr == 0)), "center")
				cov.locale <- which(chr.pos[,1] == 0)
				# ph.x[cov.locale] <- cov.centers[covar.pos]
				ph.x[cov.locale] <- cov.centers
				}
			}
		
		#sort chr.x and ph.x
		chr.x <- chr.x[order(chr.x[,1]),,drop = FALSE]
		ph.x <- sort(ph.x)
				
		# dev.new(width = 12, height = 7)
		plot.new()
		plot.window(xlim = c(0,length(chr)), ylim = c(0,1))
			
		all.pheno.y <- NULL
		par(xpd = TRUE)
		num.pheno <- length(pheno)
		plot.dim <- par("usr")
		ph.y <- (abs(plot.dim[4] - plot.dim[3]))*0.8
		ph.x.label <- plot.dim[1]
		if(num.pheno > 10){
			gap <- abs(ph.y-plot.dim[3])/num.pheno
			}else{
			gap <- abs(ph.y-plot.dim[3])/10
			}
		for(ph in 1:num.pheno){
			edge.col <- rep("gray", length(ph.x))
			names(edge.col) <- names(ph.x)
			sig.locale <- which(adj.mat[,pheno[ph]] != 0)
			if(length(sig.locale) > 0){
				pos.locale <- which(names(edge.col) %in% rownames(adj.mat)[(which(adj.mat[,pheno[ph]] > 0))])
				edge.col[pos.locale] <- rgb(0, 201/256, 87/256)
				neg.locale <- which(names(edge.col) %in% rownames(adj.mat)[(which(adj.mat[,pheno[ph]] < 0))])
				edge.col[neg.locale] <- "red"
				}
			text(x = ph.x.label, y = ph.y, pheno[ph], cex = 0.8, adj = 1)
			all.pheno.y <- c(all.pheno.y, ph.y)
			if(show.effect.size){
				segments(x0 = ph.x, y0 = rep(ph.y, length(ph.x)), y1 = rep(ph.y, length(ph.x))+(rel.block.effects[,ph]*gap*0.5), col = edge.col)	
				}else{
				points(x = ph.x, y = rep(ph.y, length(ph.x)), col = edge.col, pch = "|", cex = 1.5)
				}
			ph.y <- ph.y - gap
			}
		inter.y <- ph.y - gap #define inter.y here in case there are no interactions
		
		#replace Chr 0 with covariate names
		chr.labels <- chr
		chr.labels <- chr.labels[which(chr.labels != 0)]
		chr.cex <- rep(1.5, length(chr.labels))
		chr.srt <- rep(0, length(chr.labels))
		label.adj <- rep(0.5, length(chr.labels))
		
		chr.labels <- c(covar.names, chr.labels)
		chr.cex <- c(rep(1, length(covar.names)), chr.cex)
		chr.srt <- c(rep(90, length(covar.names)), chr.srt)
		label.adj <- c(rep(1, length(covar.names)), label.adj)

		#add the chromosome segments, labels and covariate labels above the phenotypes
		par(xpd = TRUE)
		for(i in 1:length(chr.x[,1])){
			segments(chr.x[i,1], max(all.pheno.y)*1.1, chr.x[i,2], max(all.pheno.y)*1.1, lwd = 5)
			text(x = mean(chr.x[i,]), y = max(all.pheno.y)*1.2, chr.labels[i], cex = chr.cex[i], srt = chr.srt[i], adj = label.adj[i])
			}
		par(xpd = FALSE)

		#add boxes for linkage blocks
		get.linkage.block.names <- function(block){
			marker.locale <- which(colnames(data.obj$geno) %in% block)	
			return(data.obj$marker.names[marker.locale])
			}
		linked.blocks <- data.obj$linkage.blocks.collapsed
		#restrict the linked blocks to only those in the requested chromosome
		linked.block.chr <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(names(linked.blocks), "Chr"), function(x) x[2])), "_"), function(x) x[1])))
		linked.block.chr[which(is.na(linked.block.chr))] <- 0
		chr.locale <- which(linked.block.chr %in% chr)
		
		linked.block.names <- lapply(linked.blocks, get.linkage.block.names)
		#go through each block and find the x positions, encompass all phenotypes
		for(l in chr.locale){
			marker.locale <- which(names(ph.x) %in% linked.block.names[[l]])
			x.adj <- 0.015; y.adj <- 0.04
			min.block.x <- min(ph.x[marker.locale]) - plot.dim[2]*x.adj
			max.block.x <- max(ph.x[marker.locale]) + plot.dim[2]*x.adj
			min.block.y <- min(all.pheno.y) - plot.dim[4]*y.adj
			max.block.y <- max(all.pheno.y) + plot.dim[4]*y.adj
			polygon(x = c(min.block.x, max.block.x, max.block.x, min.block.x), y = c(min.block.y, min.block.y, max.block.y, max.block.y), lwd = 2)
			}
	
		#now add colored lines for interactions
	
		plot.inter <- function(interaction.mat, y.mat, col){
			
			#replace the interacting marker locations with chromosomes
			interaction.mat[,2] <- adj.chr[(interaction.mat[,2])]
						
			#make sure all the rownames are for the chromosome in question
			rownames(interaction.mat) <- rownames(source.interactions)[as.numeric(interaction.mat[,1])]
						
			par(xpd = TRUE)
			u_inter.chr <- sort(as.numeric(unique(interaction.mat[,2])))
			marker.border <- plot.dim[2]*0.02
			line.width = 2*marker.border*interaction.line.width
			adj.marker.border <- line.width/2
			for(ch in 1:length(u_inter.chr)){
				chr.locale <- which(interaction.mat[,2] == u_inter.chr[ch])
				inter.markers <- interaction.mat[chr.locale,1,drop=FALSE]
				marker.locale <- which(names(ph.x) %in% rownames(inter.markers))
				marker.x <- ph.x[marker.locale]
				y.val <- y.mat[which(rownames(y.mat) == u_inter.chr[ch]),1]
				segments(x0 = (marker.x-adj.marker.border), y0 = rep(y.val, length(marker.locale)), x1 = (marker.x+adj.marker.border), col = col, lwd = 3)
				# text(x = plot.dim[1], y = y.val, labels = u_inter.chr[ch])
				}
			par(xpd = FALSE)
			}

		
		just.marker.locale <- dim(adj.mat)[1]
		
		#get the chromosome labels for the markers in the adjacency matrix
		adj.chr <- data.obj$chromosome[which(data.obj$marker.names %in% rownames(adj.mat))]
		
		chr.locale <- which(adj.chr %in% chr)
		source.interactions <- adj.mat[chr.locale,1:just.marker.locale]
		target.interactions <- adj.mat[1:just.marker.locale, chr.locale]
		
		#figure out which markers on the chromosome are sources of interactions, 
		#and which chromosomes they interact with. Plot the source and target
		#interactions in different colors
		any.interactions <- 0
		
		source.locale <- which(source.interactions != 0, arr.ind = TRUE)
		target.locale <- which(target.interactions != 0, arr.ind = TRUE)
		
		if((length(source.locale)+length(target.locale)) > 0){
			target.locale[,1:2] <- target.locale[,2:1] #reverse the target data to match the source data		
			inter.locale <- rbind(source.locale, target.locale)

			inter.locale[,2] <- adj.chr[(inter.locale[,2])]
			u_chr <- sort(unique(as.numeric(inter.locale[,2])))
			inter.y.mat <- matrix(NA, ncol = 1, nrow = length(u_chr))
			inter.y <- min(all.pheno.y) - gap
			inter.gap <- gap*0.5
			inter.y.mat[1,1] <- inter.y
			rownames(inter.y.mat) <- u_chr
			if(length(u_chr) > 1){
				for(i in 2:length(u_chr)){
					inter.y <- inter.y - inter.gap
					inter.y.mat[i,] <- inter.y
					}
				}
			par(xpd = TRUE)
			for(i in 1:dim(inter.y.mat)[1]){
				text(x = plot.dim[1], y = inter.y.mat[i,1], labels = rownames(inter.y.mat)[i])
				}
			par(xpd = FALSE)		
			
			if(length(source.locale) > 0){
					any.interactions <- 1
					plot.inter(source.locale, inter.y.mat, col = "blue")
					}
		
			
			inter.y.mat <- inter.y.mat-gap*0.1
			if(length(target.locale) > 0){
					any.interactions <- 1
					plot.inter(interaction.mat = target.locale, y.mat = inter.y.mat, col = "purple")
					}
			
			legend("bottomleft", legend = c("positive main effect", "negative main effect", "incoming interaction", "outgoing interaction"), lwd = 3, ncol = 2, col = c("green", "red", "purple", "blue"))
			
			}
		
		if(any.interactions == 0){
			text(mean(plot.dim[1:2]), y = inter.y, labels = "No Interactions")
			}
		}
	}