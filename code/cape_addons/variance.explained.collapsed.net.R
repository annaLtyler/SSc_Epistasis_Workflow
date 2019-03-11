#This is a new variance explained 
#color.scheme can either be "grays"
#or the name of a ColorBrewer pallete
#data.obj <- readRDS("~/Documents/Data/DO/Results/all_pheno/crossDO.RData")

variance.explained.collapsed.net <- function(data.obj, color.scheme = c("Paired", "grays")){

	library(RColorBrewer)	
	color.scheme <- color.scheme[1]
	
	link.blocks <- data.obj$linkage.blocks.collapsed
#====================================================================
#internal functions
#====================================================================
get.max.main <- function(block.name, pheno.name){
	block.markers <- link.blocks[[which(names(link.blocks) == block.name)]]
	block.locale <- which(colnames(data.obj$geno.for.pairscan) %in% block.markers)
	block.effects <- data.obj$full.net[block.locale, pheno.name,drop=FALSE]
	max.block <- which.max(abs(block.effects))
	return(data.obj$linkage.blocks.full[[rownames(block.effects)[max.block]]])
	}

get.max.int <- function(source.block, target.block){
	source.markers <- link.blocks[[which(names(link.blocks) == source.block)]]
	source.locale <- which(colnames(data.obj$geno.for.pairscan) %in% source.markers)
	target.markers <- link.blocks[[which(names(link.blocks) == target.block)]]
	target.locale <- which(colnames(data.obj$geno.for.pairscan) %in% target.markers)
	int.mat <- data.obj$full.net[source.locale, target.locale,drop=FALSE]
	max.locale <- which(abs(int.mat) == max(abs(int.mat)), arr.ind = TRUE)
	max.source <- data.obj$linkage.blocks.full[[rownames(int.mat)[max.locale[1,1]]]]
	max.target <- data.obj$linkage.blocks.full[[colnames(int.mat)[max.locale[1,2]]]]
	return(c(max.source, max.target))
	}
#====================================================================
	
	
	covar.info <- get.covar(data.obj)
	
	
	net <- data.obj$collapsed.net
	just.int <- net[,1:nrow(net)]
	just.main <- net[,(nrow(net)+1):ncol(net)]
	
	if(length(covar.info$covar.names) > 0){
		covar.locale <- match(covar.info$covar.names, rownames(just.int))
		just.int <- just.int[-covar.locale,-covar.locale]
		just.main <- just.main[-covar.locale,]
		covar.table <- covar.info$covar.table
		colnames(covar.table) <- covar.info$covar.names
		}else{
		covar.table <- NULL
		}

	all.marker.names <- colnames(data.obj$geno.for.pairscan)
	
	#create unique names for each marker that doesn't have any numbers in it
	#the number names can't be used in a formula
	dummy.marker.names <- sapply(1:length(all.marker.names), function(x) paste(sample(c(letters, LETTERS), 5), collapse = ""))
	
	
	pheno <- get.pheno(data.obj, "raw.traits", covar = NULL)

	bare.models <- matrix(NA, nrow = length(covar.info$covar.names), ncol = ncol(pheno))
	rownames(bare.models) <- covar.info$covar.names
	
	for(cv in 1:length(covar.info$covar.names)){
		cv.model <- apply(pheno, 2, function(x) lm(x~data.obj$p.covar.table[,cv]))		
		bare.models[cv,] <- sapply(cv.model, function(x) summary(x)$adj.r.squared)
		}
	
	
	main.r2 <- rep(NA, ncol(pheno))
	full.r2 <- rep(NA, ncol(pheno))
	num.main <- rep(NA, ncol(pheno))
	num.int <- rep(NA, ncol(pheno))
	names(main.r2) <- names(full.r2) <- names(num.main) <- names(num.int) <- colnames(pheno)
	
	geno.data <- data.obj$geno.for.pairscan
	colnames(geno.data) <- dummy.marker.names
	
	data.fm <- data.frame(cbind(pheno, covar.table, geno.data))
	
	for(ph in 1:ncol(pheno)){
		#find each marker with a main effect on the phenotype
		main.locale <- which(just.main[,ph] != 0)
		num.main[ph] <- length(main.locale)
		
		max.marker.names <- sapply(names(main.locale), function(x) get.max.main(x, colnames(pheno)[ph]))
		max.marker.locale <- which(colnames(data.obj$geno.for.pairscan) %in% max.marker.names)
		
		#make a model with all of the main effects
		marker.names <- dummy.marker.names[max.marker.locale]

		#because we look at the collapsed network, we need to 
		#find the individual markers that best represent each
		#block

		main.geno <- geno.data[,max.marker.locale]
		
		main.model <- lm(pheno[,ph]~data.obj$p.covar.table+main.geno)
		main.r2[ph] <- summary(main.model)$adj.r.squared

		#now find all interactions that each main effect marker takes
		#part in and build a model with all main effects and interactions
		covar.text <- paste(covar.info$covar.names, collapse = "+")
		marker.text <- paste(marker.names, collapse = "+")
		fmla <- paste(colnames(pheno)[ph], "~", covar.text, "+", marker.text)
		
		int.table <- NULL
		
		for(m in 1:length(main.locale)){
			
			#get all interaction values for the source markers targeting
			#the main effect marker
			source.locale <- which(just.int[,main.locale[m]] != 0)
			if(length(source.locale) > 0){
				for(sc in 1:length(source.locale)){
					source.block <- names(source.locale)[sc]
					target.block <- names(main.locale)[m]
					max.markers <- get.max.int(source.block, target.block)
					max.marker.locale <- match(max.markers, colnames(data.obj$geno.for.pairscan))
					int.table <- rbind(int.table, c(dummy.marker.names[max.marker.locale[1]], dummy.marker.names[max.marker.locale[2]]))
					}
				}

			target.locale <- which(just.int[main.locale[m],] != 0)
			if(length(target.locale) > 0){
				for(tg in 1:length(target.locale)){
					target.block <- names(target.locale)[tg]
					source.block <- names(main.locale)[m]
					max.markers <- get.max.int(source.block, target.block)
					max.marker.locale <- match(max.markers, colnames(data.obj$geno.for.pairscan))
					int.table <- rbind(int.table, c(dummy.marker.names[max.marker.locale[1]], dummy.marker.names[max.marker.locale[2]]))
					}

				}

			} #end looping through all main effect markers for this phenotype

		#sort the table so we only include unique interactions
		if(length(int.table) > 0){
			sorted.table <- unique(t(apply(int.table, 1, sort)))
			int.text <- paste(apply(sorted.table, 1, function(x) paste(x[1], x[2], sep = ":")), collapse = "+")
			final.fmla <- paste(fmla, "+", int.text)
			num.int[ph] <- nrow(int.table)
			}else{
			final.fmla <- fmla
			num.int[ph] <- 0
			}
		

		full.model <- lm(as.formula(final.fmla), data = data.fm)
		full.r2[ph] <- summary(full.model)$adj.r.squared

		}
	
	final.results.table <- rbind(bare.models, main.r2, num.main, full.r2, num.int)
	final.results.table[which(final.results.table < 0)] <- 0
	bare.models[which(bare.models < 0)] <- 0
	
	results.to.plot <- rbind(bare.models, main.r2, full.r2)
	covar.text <- paste(covar.info$covar.names, collapse = " + ")
	# cat("plotting to Variance.Explained.Collapsed.Net.pdf\n")
	
	# pdf("Variance.Explained.Collapsed.Net.pdf", width = 10, height = 6)
	
	if(color.scheme == "grays"){
	cols <- gray.colors(nrow(results.to.plot))
		}else{
		cols <- brewer.pal(nrow(results.to.plot), color.scheme)
		}
		
	# par(mar = c(6, 4, 2, 2))
	# barplot(results.to.plot, beside = TRUE, col = cols[1:nrow(results.to.plot)], ylim = c(0, 1), main = "Variance Explained")
	# legend("topleft", fill = cols[1:nrow(results.to.plot)], legend = c(covar.info$covar.names, "main effects", "interactions"))
	# plot.dim <- par("usr")
	# plot.height <- plot.dim[4] - plot.dim[3]
	# plot.width <- plot.dim[2] - plot.dim[1]
	
	# par(xpd = TRUE)
	# text(x = segment.region(plot.dim[1], plot.dim[2], ncol(results.to.plot)), y = rep(plot.dim[3]-(plot.height*0.15), ncol(results.to.plot)), labels = num.main)
	# text(x = (plot.dim[1] - (plot.width*0.01)), y = (plot.dim[3]-(plot.height*0.15)), labels = "num Main")
	
	# text(x = segment.region(plot.dim[1], plot.dim[2], ncol(results.to.plot)), y = rep(plot.dim[3]-(plot.height*0.2), ncol(results.to.plot)), labels = num.int)
	# text(x = (plot.dim[1] - (plot.width*0.01)), y = (plot.dim[3]-(plot.height*0.2)), labels = "num Int")
	# par(xpd = FALSE)
	

	gained.by.main <- main.r2 - colSums(bare.models)
	gained.by.int <- full.r2 - main.r2
	
	amount.gained <- rbind(bare.models, gained.by.main, gained.by.int)
	amount.gained[which(amount.gained < 0)] <- 0
	max.ex <- max(colSums(amount.gained))
	#round up to nearest 10th
	maxx <- ceiling(max.ex*10)/10

	layout(matrix(c(1,2), ncol = 1), heights = c(0.2, 1))
	par(mar = c(0,0,0,0))
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	legend(0,1, fill = cols[1:nrow(results.to.plot)], legend = c(covar.info$covar.names, "Main Effects", "Interactions"), cex = 0.9)
	par(mar = c(4,6,3,4))
	barplot(amount.gained, beside = FALSE, main = "Variance Explained by Main Effects and Interactions", col = cols[1:nrow(results.to.plot)], horiz = TRUE, las = 2, adj = 0, xlim = c(0, maxx))
	# dev.off()
	
	# final.results.table <- rbind(final.results.table, amount.gained)
	return(amount.gained)


	}