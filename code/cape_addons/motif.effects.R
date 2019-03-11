#This function looks at the main effects in different motifs

motif.effects <- function(data.obj, interaction = c("enhancing", "suppressing"), main = c("coherent", "incoherent"), source.sign = c("pos", "neg", "any")){
	
	# data.obj <- get.network(data.obj, p.or.q = p.or.q, collapse.linked.markers = TRUE)
	motif.obj <- find.motifs(data.obj)
	
	#collect the unique motifs for all phenotypes
	motifs <- NULL
	for(i in 1:length(motif.obj[[2]])){
		if(interaction == "enhancing"){
			motif.locale.int <- which(motif.obj[[2]][[i]][,1] == 1)
			}else{
			motif.locale.int <- which(motif.obj[[2]][[i]][,1] == -1)	
			}
		if(main == "coherent"){
			motif.locale.main <- which(motif.obj[[2]][[i]][,2] == motif.obj[[2]][[i]][,3])
			}else{
			motif.locale.main <- which(motif.obj[[2]][[i]][,2] != motif.obj[[2]][[i]][,3])	
			}
		if(source.sign == "pos"){
			sign.locale <- which(motif.obj[[2]][[i]][,2] == 1)
			}
		if(source.sign == "neg"){
			sign.locale <- which(motif.obj[[2]][[i]][,2] == -1)
			}
		if(source.sign == "any"){
			sign.locale <- 1:dim(motif.obj[[2]][[i]])[1]
			}
		motif.locale <- intersect(intersect(motif.locale.int, motif.locale.main), sign.locale)
		motifs <- rbind(motifs, motif.obj[[1]][[i]][motif.locale,1:2])
		}
	
	u_motifs <- unique(motifs)
	
	#now fill two matrices. one with the main effects of the source
	#marker and one with the main effects of the target marker
	
	collapsed.net <- data.obj$collapsed.net
	num.markers <- min(dim(collapsed.net))
	just.pheno <- collapsed.net[,(num.markers+1):dim(collapsed.net)[2]]
	
	source.locale <- match(u_motifs[,1], rownames(collapsed.net))
	source.mat <- just.pheno[source.locale,]
	source.effect.mat <- matrix(0, nrow = nrow(source.mat), ncol = ncol(source.mat))
	source.effect.mat[which(source.mat != 0)] <- 1
	source.sum <- rowSums(source.effect.mat)
	
	target.locale <- match(u_motifs[,2], rownames(collapsed.net))
	target.mat <- just.pheno[target.locale,]
	target.effect.mat <- matrix(0, nrow = nrow(target.mat), ncol = ncol(target.mat))
	target.effect.mat[which(target.mat != 0)] <- 1
	target.sum <- rowSums(target.effect.mat)
	
	#cluster by the target phenotype effects, since this is a large part of how directionality is determined
	tm <- heatmap(target.mat)
	
	plot.mat <- function(mat){
		ymin <- min(mat);ymax <- max(mat)
		xmin <- 1;xmax <- ncol(mat)
		plot.height <- ymax - ymin
		par(mar = c(2, 2, 2, 2))
		plot.new()
		plot.window(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
		for(i in 1:ncol(mat)){
			points(x = jitter(rep(i, nrow(mat))), mat[,i])
			}
		axis(2)
		par(xpd = TRUE)
		text(x = 1:ncol(mat), y = rep(ymin-plot.height*0.1, ncol(mat)), labels = colnames(mat))
		par(xpd = FALSE)
		}
	
	# quartz()
	layout.mat <- matrix(c(1:4), ncol = 2)
	layout(layout.mat, heights = c(1,0.5))
	par(mar = c(4,2,3,2))
	
	imageWithText(source.mat[tm$rowInd,], show.text = FALSE, split.at.vals = TRUE, col.names = names(motif.obj[[1]]), col.text.adj = 1, col.text.shift = 0, main = paste("Source\n", interaction, main, "\nSource =", source.sign))
	par(mar = c(3,2,0,2))
	plot.mat(source.mat)
	par(mar = c(4,2,3,2))
	imageWithText(target.mat[tm$rowInd,], show.text = FALSE, split.at.vals = TRUE, col.names = names(motif.obj[[1]]), col.text.adj = 1, col.text.shift = 0, main = paste("Target\n", interaction, main, "\nSource =", source.sign))
	par(mar = c(3,2,0,2))
	# boxplot(target.mat, las = 2)
	plot.mat(target.mat)

	ymin <- min(c(source.mat, target.mat))
	ymax <- max(c(source.mat, target.mat))

	
	# layout.mat <- matrix(c(1:4), ncol = 2)
	# layout(layout.mat)
	# for(ph in 1:dim(source.mat)[2]){
		# plot.new()
		# plot.window(xlim = c(0,1), ylim = c(ymin, ymax))
		# segments(x0 = rep(0.2, dim(source.mat)[1]), y0 = source.mat[,ph], x1 = rep(0.8, dim(source.mat)[1]), y1 = target.mat[,ph])
		# segments(0,0,1,0, col = "red")
		# mtext(names(motif.obj[[2]])[ph], cex = 2)
		# text(x = 0.2, y = ymin-0.15, "Source")
		# text(x = 0.8, y = ymin-0.15, "Target")
		# mtext(paste(interaction, main, "\nSource =", source.sign), outer = TRUE, line = -2)
		# }

	# quartz()
	# imageWithText(cbind(source.sum, target.sum, (source.sum-target.sum)), show.text = FALSE, split.at.vals = TRUE)	
	
	
}