#This is a new variance explained 
#color.scheme can either be "grays"
#or the name of a ColorBrewer pallete


variance.explained <- function(data.obj, color.scheme = c("grays", "Paired")){

	library(RColorBrewer)	
	if(length(color.scheme) > 1){
		color.scheme = "grays"
		}
	
	covar.info <- get.covar(data.obj)
	
	
	net <- data.obj$full.net
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
		
		#make a model with all of the main effects
		marker.names <- dummy.marker.names[main.locale]

		main.geno <- geno.data[,main.locale]
		
		main.model <- lm(pheno[,ph]~data.obj$p.covar.table+main.geno)
		main.r2[ph] <- summary(main.model)$adj.r.squared

		#now find all interactions that each main effect marker takes
		#part in and build a model with all main effects and interactions
		covar.text <- paste(covar.info$covar.names, collapse = "+")
		marker.text <- paste(marker.names, collapse = "+")
		fmla <- paste(colnames(pheno)[ph], "~", covar.text, "+", marker.text)
		
		int.table <- NULL
		
		for(m in 1:length(main.locale)){
			
			source.locale <- which(just.int[,main.locale[m]] != 0)
			if(length(source.locale) > 0){
				for(sc in 1:length(source.locale)){
					int.table <- rbind(int.table, c(dummy.marker.names[source.locale[sc]], dummy.marker.names[main.locale[m]]))
					}
				}

			target.locale <- which(just.int[main.locale[m],] != 0)
			if(length(target.locale) > 0){
				for(tg in 1:length(target.locale)){
					int.table <- rbind(int.table, c(dummy.marker.names[main.locale[m]], dummy.marker.names[target.locale[tg]]))
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
	cat("plotting to Variance.Explained.pdf\n")
	
	pdf("Variance.Explained.pdf", width = 10, height = 6)
	
	if(color.scheme == "grays"){
	cols <- gray.colors(nrow(results.to.plot))
		}else{
		cols <- brewer.pal(nrow(results.to.plot), color.scheme)
		}
		
	par(mar = c(6, 4, 2, 2))
	barplot(results.to.plot, beside = TRUE, col = cols[1:nrow(results.to.plot)], ylim = c(0, 1), main = "Variance Explained")
	legend("topleft", fill = cols[1:nrow(results.to.plot)], legend = c(covar.info$covar.names, "main effects", "interactions"))
	plot.dim <- par("usr")
	plot.height <- plot.dim[4] - plot.dim[3]
	plot.width <- plot.dim[2] - plot.dim[1]
	
	par(xpd = TRUE)
	text(x = segment.region(plot.dim[1], plot.dim[2], ncol(results.to.plot)), y = rep(plot.dim[3]-(plot.height*0.15), ncol(results.to.plot)), labels = num.main)
	text(x = (plot.dim[1] - (plot.width*0.01)), y = (plot.dim[3]-(plot.height*0.15)), labels = "num Main")
	
	text(x = segment.region(plot.dim[1], plot.dim[2], ncol(results.to.plot)), y = rep(plot.dim[3]-(plot.height*0.2), ncol(results.to.plot)), labels = num.int)
	text(x = (plot.dim[1] - (plot.width*0.01)), y = (plot.dim[3]-(plot.height*0.2)), labels = "num Int")
	par(xpd = FALSE)
	

	gained.by.main <- main.r2 - colSums(bare.models)
	gained.by.int <- full.r2 - main.r2
	
	amount.gained <- rbind(bare.models, gained.by.main, gained.by.int)
	max.ex <- max(colSums(amount.gained))
	#round up to nearest 10th
	maxx <- ceiling(max.ex*10)/10

	layout(matrix(c(1,2), ncol = 1), heights = c(0.2, 1))
	par(mar = c(0,0,0,0))
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	legend(0,1, fill = cols[1:nrow(results.to.plot)], legend = c(covar.info$covar.names, "main effects", "interactions"), cex = 0.9)
	par(mar = c(4,5,3,4))
	barplot(amount.gained, beside = FALSE, main = "Variance Explained by Main Effects and Interactions", col = cols[1:nrow(results.to.plot)], horiz = TRUE, las = 2, adj = 0, xlim = c(0, maxx))
	dev.off()
	
	return(final.results.table)


	}