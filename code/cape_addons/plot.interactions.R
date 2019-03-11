#This function plots significant interactions


plot.interactions <- function(data.obj, chr = NULL, Lblock = NULL, p.or.q = 0.05, pdf.label = "Interaction.Plots.pdf", error.bars = TRUE, normalized.pheno = TRUE, bin.geno = c(0,0.5,1), covar = NULL, alt.bin.var = NULL, alt.bin.bins = NULL){
	
	var.influences <- data.obj$var.to.var.p.val
	
	get.chr <- function(marker.name){
		marker.locale <- which(colnames(data.obj$geno) == marker.name)
		return(data.obj$chromosome[marker.locale])
		}

	get.block <- function(marker.name){
		return(names(unlist(lapply(data.obj$linkage.blocks.collapsed, function(x) which(x == marker.name)))))	
		}

	if(!is.null(chr)){
		cat("Finding marker chromosomes...\n")
		marker.chr <- t(apply(var.influences[,1:2], 1, function(x) c(get.chr(x[1]), get.chr(x[2]))))
		}
	if(!is.null(Lblock)){
		cat("Finding marker linkage blocks...\n")
		marker.block <- t(apply(var.influences[,1:2], 1, function(x) c(get.block(x[1]), get.block(x[2]))))
		}
	
	pheno.names <- colnames(data.obj$pheno)
	num.pheno <- length(pheno.names)
	num.et <- length(data.obj$ET[1,])
	et.names <- colnames(data.obj$ET)
	geno <- data.obj$geno
	
	
	if(normalized.pheno){
		all.pheno <- cbind(data.obj$pheno, data.obj$ET)
		}else{
		#use the uncentered, unnormalized data
		all.pheno <- cbind(data.obj$raw.pheno, data.obj$ET)
		}
		
	all.pheno.names <- c(pheno.names, et.names)
	total.response <- num.et + num.pheno
	
	if(is.null(var.influences)){
		stop("variant-to-variant influences must be calculated before plotting.")
		}


	chr.locale <- NULL; block.locale <- NULL
	#subset out only the chromosomes or blocks asked for
	if(!is.null(chr)){
		chr.locale <- unique(c(which(marker.chr[,1] %in% chr), which(marker.chr[,2] %in% chr)))
		}
	if(!is.null(chr)){
		block.locale <- unique(c(which(marker.block[,1] %in% Lblock), which(marker.block[,2] %in% Lblock)))
		}

	subset.locale <- unique(c(chr.locale, block.locale))
	if(length(subset.locale) > 0){
		var.influences <- var.influences[subset.locale,]
		}

	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))
	sig.int.locale <- which(var.influences[,var.sig.col] <= p.or.q)
	
	if(length(sig.int.locale) == 0){
		
		pdf(pdf.label)
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(0,1))
		text(0.5, 0.5, "No Significant Interactions")
		dev.off()
		}else{
		
		#take out the covariate from the interaction list
		if(!is.null(covar)){
			marker1.names <- data.obj$marker.names[match(var.influences[,1], colnames(data.obj$geno))]
			marker2.names <- data.obj$marker.names[match(var.influences[,2], colnames(data.obj$geno))]
			covar.locale <- intersect(c(which(marker1.names %in% covar), which(marker2.names %in% covar)), sig.int.locale)
			if(length(covar.locale) > 0){
				sig.int.locale <- sig.int.locale[-match(covar.locale, sig.int.locale)]
				}			
			}
		
		num.int <- length(sig.int.locale)
		if(num.int > 4){
			num.per.page <- 4
			}else{
				num.per.page <- num.int
				}
	
		plot.obj.list <- vector(length = length(sig.int.locale), mode = "list")
		plot.obj.count <- 1
		cat("calculating interactions...\n")
			for(i in 1:length(sig.int.locale)){
				report.progress(i, length(sig.int.locale))
				plot.obj <- vector(length = 8, mode = "list")
				names(plot.obj) <- c("marker1", "marker2", "pheno", "main", "xlab", "ylab", "trace.label", "errors")
			
				marker1.name <- data.obj$marker.names[which(colnames(geno) == var.influences[sig.int.locale[i],1])]
				marker2.name <- data.obj$marker.names[which(colnames(geno) == var.influences[sig.int.locale[i],2])]
				marker1.block <- names(unlist(lapply(data.obj$linkage.blocks.collapsed, function(x) which(x ==  var.influences[sig.int.locale[i],1]))))
				marker2.block <- names(unlist(lapply(data.obj$linkage.blocks.collapsed, function(x) which(x ==  var.influences[sig.int.locale[i],2]))))
				inter.p <- signif(var.influences[sig.int.locale[i],7], 3)
				inter.effect <- signif(var.influences[sig.int.locale[i],3]/var.influences[sig.int.locale[i],4], 3)
				
				if(!(marker1.name %in% covar)){
					plot.obj$marker1 <- bin.vector(geno[,as.character(var.influences[sig.int.locale[i],1])], bin.geno)
					}
				if(!(marker2.name %in% covar)){
					plot.obj$marker2 <- bin.vector(geno[,as.character(var.influences[sig.int.locale[i],2])], bin.geno)	
					}

				if((marker1.name %in% alt.bin.var)){
					plot.obj$marker1 <- bin.vector(geno[,as.character(var.influences[sig.int.locale[i],1])], alt.bin.bins)	
					}
				if((marker2.name %in% alt.bin.var)){
					plot.obj$marker2 <- bin.vector(geno[,as.character(var.influences[sig.int.locale[i],2])], alt.bin.bins)	
					}

				plot.obj$xlab = marker1.name
				plot.obj$trace.label = marker2.name

				for(ph in 1:total.response){
					plot.obj$ylab = all.pheno.names[ph]
					#if we are adjusting for a covariate, do this here
					if(!is.null(covar)){
						covar.col <- match(covar, data.obj$marker.names)
						covar.v <- data.obj$geno[,covar.col,drop = FALSE]
						covar.means <- apply(covar.v, 2, function(x) mean(x, na.rm = TRUE))
						for(i in 1:length(covar)){
							covar.v[,i] <- covar.v[,i] - covar.means[i]
							}
						model <- lm(all.pheno[,ph]~covar.v)
						plot.obj$pheno <- residuals(model)
						
						na.locale <- unique(which(is.na(cbind(all.pheno[,ph], geno[,covar.col])), arr.ind = TRUE)[,1])
						not.na <- c(1:dim(all.pheno)[1])
						if(length(na.locale) > 0){
							not.na <- not.na[-na.locale]
							}
						plot.obj$not.na <- not.na
						
						}else{
						plot.obj$pheno <- all.pheno[,ph]
						na.locale <- unique(which(is.na(all.pheno[,ph])))
						not.na <- c(1:dim(all.pheno)[1])
						if(length(na.locale) > 0){
							not.na <- not.na[-na.locale]
							}
						plot.obj$not.na <- not.na

						}
					plot.obj$main <- paste(marker1.block, ":", marker1.name, ", ", marker2.block, ":", marker2.name, "\n", all.pheno.names[ph], " Eff:", inter.effect, " p:", inter.p, sep = "")

					plot.obj$errors <- get.interaction.error(plot.obj$marker1[not.na], plot.obj$marker2[not.na], plot.obj$pheno, error.type = "sd")
	
					plot.obj.list[[plot.obj.count]] <- plot.obj
					plot.obj.count <- plot.obj.count + 1
					}					
				}
				
				
			get.pheno.ylim <- function(plot.obj){
				pheno.name <- plot.obj$ylab
				if(error.bars){
					all.ylim <- c(min(plot.obj$errors$means - plot.obj$errors[[2]], na.rm = TRUE), max(plot.obj$errors$means + plot.obj$errors[[2]], na.rm = TRUE))
					}else{
					all.ylim <- c(min(plot.obj$errors$means, na.rm = TRUE), max(plot.obj$errors$means, na.rm = TRUE))
					}
					result <- c(pheno.name, all.ylim)
					return(result)
				}
			
			all.ylim <- lapply(plot.obj.list, get.pheno.ylim)		
			
			overall.ylim <- matrix(NA, nrow = length(all.pheno.names), ncol = 2)
			for(p in 1:length(all.pheno.names)){
				all.pheno.lim <- unlist(lapply(all.ylim, function(x) if(x[1] == all.pheno.names[p]){return(x[2:3])}))
				overall.ylim[p,] <- c(min(as.numeric(all.pheno.lim)), max(as.numeric(all.pheno.lim)))
				}
				
			
			plot.interaction <- function(plot.obj){
				pheno.locale <- which(all.pheno.names %in% plot.obj$ylab)
				not.na <- plot.obj$not.na
				interaction.plot(plot.obj$marker1[not.na], plot.obj$marker2[not.na], plot.obj$pheno, ylab = plot.obj$ylab, main = plot.obj$main, xlab = plot.obj$xlab, trace.label = plot.obj$trace.label, ylim = overall.ylim[pheno.locale,], fixed = TRUE)
				if(error.bars){
					for(j in 1:length(plot.obj$errors$means[,1])){
						segments(1:length(colnames(plot.obj$errors$means)), (plot.obj$errors$means[j,]+plot.obj$errors[[2]][j,]), 1:length(colnames(plot.obj$errors$means)), (plot.obj$errors$means[j,]-plot.obj$errors[[2]][j,]))
						}
					}
				}
				
				
			cat("\nplotting results...\n")
			pdf(pdf.label, width = total.response*3, height = num.per.page*3)
			layout(matrix(1:(total.response*num.per.page), ncol = total.response, byrow = TRUE))
			# for(i in 1:length(plot.obj.list)){
				# print(i)
				# plot.interaction(plot.obj.list[[i]])	
				# }
			null <- lapply(plot.obj.list, plot.interaction)
			dev.off()

						
			}

}

