	#This code will plot null and observed distributions for
	#a data object called "data.obj"
	#It plots distributions for both the standardized interaction
	#coefficients and m12/m21
	#if you plot a main effect null, we actually have to
	#calculate it, so it will take a long time

plot.null.dist <- function(data.obj, pairscan.obj, path = "."){
	require(timeDate)
	require(GGally)
	
	if(missing(pairscan.obj)){
		pairscan.results <- data.obj$pairscan.results	
		pairscan.perm <- data.obj$pairscan.perm		
		if(is.null(pairscan.results)){
			stop("I cannot find the pairscan results in the data object.")
			}
		}else{
		pairscan.results <- pairscan.obj$pairscan.results	
		pairscan.perm <- pairscan.obj$pairscan.perm
		}
		if(is.null(pairscan.results)){
			stop("I cannot find the pairscan results.")
			}		
	#====================================================================================================
	# interneal functions
	#====================================================================================================
	plot.dist <- function(null.dist = NULL, true.dist = NULL, plot.title = NULL, xlim = NULL, legend.cex = 2, legend.pos = "topright", offset = 0){
		if(!is.null(null.dist)){
			null.dens <- density(null.dist[which(!is.na(null.dist))])
			}else{
			null.dens <- NULL	
			}
		if(!is.null(true.dist)){
			true.dens <- density(true.dist[which(!is.na(true.dist))])
			}else{
			true.dens <- NULL	
			}
		if(is.null(xlim)){
			x.min <- min(null.dens$x, true.dens$x); x.max <- max(null.dens$x, true.dens$x)
			xlim <- c(x.min, x.max)
			}else{
			x.min <- xlim[1]; x.max <- xlim[2]	
			}
		cols <- c(rgb(190/256,174/256,212/256), rgb(127/256,201/256,127/256))
		y.min <- min(null.dens$y, true.dens$y); y.max <- max(null.dens$y+offset, true.dens$y+offset)
		plot(null.dens, main = plot.title, col = cols[1], xlim = c(x.min, x.max), ylim = c(y.min, y.max), lwd = 3, cex.lab = 1.5, cex.axis = 2)
		points(x = true.dens$x, y = true.dens$y+offset, col = cols[2], type = "l", lwd = 3)
		abline(v = median(null.dist), lwd = 1)
		legend(legend.pos, legend = c("Null", "Obs."), fill = cols, cex = legend.cex)		
		}
			
		
	#distribution of m12/m21
	# layout.mat <- matrix(1:12, ncol = 4, byrow = TRUE)
	layout.mat <- matrix(c(1:8,0), ncol = 3, byrow = TRUE)
	jpeg(paste0(path, "/null.distributions.m.jpg"), width = dim(layout.mat)[2]*400, height = dim(layout.mat)[1]*400, quality = 100)
	layout(layout.mat)
	# layout.show(8)
	
	
	#================================================================================
	# Plot the -log10 of the m distributions, null and observed
	#================================================================================
	null.m12 <- log10(abs(as.numeric(data.obj$var.to.var.influences.perm[,"m12"])))
	null.m21 <- log10(abs(as.numeric(data.obj$var.to.var.influences.perm[,"m21"])))
	
	true.m12 <- log10(abs(as.numeric(data.obj$var.to.var.influences[,"m12"])))
	true.m21 <- log10(abs(as.numeric(data.obj$var.to.var.influences[,"m21"])))
	
	null.m <- c(null.m12, null.m21)
	true.m <- c(true.m12, true.m21)
	
	plot.dist(null.m, true.m, plot.title = paste("Distribution of Log10 Raw |m|\nMedian:", signif(median(null.m), 3), "\nnote: observed distribution is offset"), legend.cex = 1.5, offset = 0.05)
	 
	#================================================================================
	# Plot the -log10 of the standard deviation of the m's
	#================================================================================	 
	null.error.m12 <- log10(as.numeric(data.obj$var.to.var.influences.perm[,"m12.std.dev"]))
	null.error.m21 <- log10(as.numeric(data.obj$var.to.var.influences.perm[,"m21.std.dev"]))
	null.error <- c(null.error.m12, null.error.m21)
	
	true.error.m12 <- log10(as.numeric(data.obj$var.to.var.influences[,"m12.std.dev"]))
	true.error.m21 <- log10(as.numeric(data.obj$var.to.var.influences[,"m21.std.dev"]))
	true.error <- c(true.error.m12, true.error.m21)

	plot.dist(null.error, true.error, plot.title = paste("Distribution of log10 m Errors:", signif(median(null.error), 3), "\nnote: observed distribution is offset"), offset = 0.05)


	#================================================================================
	# Plot the -log of the standardized m's, true and null
	#================================================================================
		
	null.std.m12 <- as.numeric(data.obj$var.to.var.influences.perm[,"m12"])/as.numeric(data.obj$var.to.var.influences.perm[,"m12.std.dev"])
	null.std.m21 <- as.numeric(data.obj$var.to.var.influences.perm[,"m21"])/as.numeric(data.obj$var.to.var.influences.perm[,"m21.std.dev"])
	null.std <- c(null.std.m12, null.std.m21)	

	true.std.m12 <- as.numeric(data.obj$var.to.var.influences[,"m12"])/as.numeric(data.obj$var.to.var.influences[,"m12.std.dev"])
	true.std.m21 <- as.numeric(data.obj$var.to.var.influences[,"m21"])/as.numeric(data.obj$var.to.var.influences[,"m21.std.dev"])
	true.std <- c(true.std.m12, true.std.m21)	
	
	plot.dist(null.std, true.std, plot.title = paste("Distribution of Standardized m\nMedian:", signif(median(null.m), 3), "\nnote: observed distribution is offset"), legend.pos = "topleft", offset = 0.2)
	
	
	#================================================================================
	#color the points in the observed distribution based on significance
	#================================================================================	
	low.null <- null.std[which(null.std < 0)]; low.fun <- ecdf(abs(low.null))
	high.null <- null.std[which(null.std > 0)]; high.fun <- ecdf(high.null)

	all.emp.p <- matrix(NA, ncol = 1, nrow = length(true.std))
	high.m.locale <- which(true.std > 0)
	low.m.locale <- which(true.std < 0)

	low.p <- 1-low.fun(abs(true.std[low.m.locale]))
	high.p <- 1-high.fun(abs(true.std[high.m.locale]))

	all.p <- rep(NA, length(true.std))
	all.p[low.m.locale] <- low.p
	all.p[high.m.locale] <- high.p

	p.col <- rep(NA, length(true.std))
	low.p.col <- colors.from.values(-log10(low.p+1e-6), grad.dir = "high", col.scale = "purple", light.dark = "d")	
	high.p.col <- colors.from.values(-log10(high.p+1e-6), grad.dir = "high", col.scale = "green", light.dark = "d")

	p.col[low.m.locale] <- low.p.col
	p.col[high.m.locale] <- high.p.col
		
	plot(true.m[low.m.locale], true.error[low.m.locale], xlab = "log10 m", ylab = "log10 m error", main = "log10 |Negative m| vs. log10 Error\nColored by -log10 p value", col = p.col[low.m.locale], pch = 16, cex.axis = 2, cex.lab = 2)
	
	
	plot(true.m[high.m.locale], true.error[high.m.locale], xlab = "log10 m", ylab = "log10 m error", main = "log10 |Positive m| vs. log10 Error\nColored by -log10 p value", col = p.col[high.m.locale], pch = 16, cex.axis = 2, cex.lab = 2)

	plot(true.m, true.std, col = p.col, pch = 16, xlab = "log10 |m|", ylab = "Std m", main = "log10 |m| vs. Std |m|", cex.axis = 2, cex.lab = 2)

	plot(true.m[low.m.locale], true.error[low.m.locale], xlab = "log10 m", ylab = "log10 m error", main = "Zoomed log10 |Negative m| vs. log10 Error\nColored by -log10 p value", col = p.col[low.m.locale], pch = 16, xlim = c(-1,1), ylim = c(-1,1), cex.axis = 2, cex.lab = 2)


	plot(true.m[high.m.locale], true.error[high.m.locale], xlab = "log10 m", ylab = "log10 m error", main = "Zoomed log10 |Positive m| vs. log10 Error\nColored by -log10 p value", col = p.col[high.m.locale], pch = 16, xlim = c(-1,1), ylim = c(-1,1), cex.axis = 2, cex.lab = 2)
	
	dev.off()
	
	
	#distribution of pairscan main effect coefficients
	layout.mat <- get.layout.mat(length(pairscan.perm))
	marker1.pos <- dim(pairscan.perm[[1]][[1]])[2] - 2
	marker2.pos <- dim(pairscan.perm[[1]][[1]])[2] - 1
	jpeg(paste0(path, "/null.distributions.main.effects.jpg"), width = dim(layout.mat)[2]*400, height = dim(layout.mat)[1]*400, quality = 100)
	layout(layout.mat)
	for(i in 1:length(pairscan.perm)){
		null.dist <- c(as.numeric(pairscan.perm[[i]][[1]][,marker1.pos]), as.numeric(pairscan.perm[[i]][[1]][,marker2.pos]))/c(as.numeric(pairscan.perm[[i]][[2]][,marker1.pos]), as.numeric(pairscan.perm[[i]][[2]][,marker2.pos]))
		true.dist <- c(as.numeric(pairscan.results[[i]][[1]][,marker1.pos]), as.numeric(pairscan.results[[i]][[1]][,marker2.pos]))/c(as.numeric(pairscan.results[[i]][[2]][,marker1.pos]), as.numeric(pairscan.results[[i]][[2]][,marker2.pos]))
		plot.dist(null.dist = abs(null.dist), true.dist = abs(true.dist), plot.title = paste("Null Distribution of\n|Standardized Main Effects| for", names(pairscan.results)[i], "\nnote observed distribution is offset"), offset = 0.05)
		}
	dev.off()
	
}