#This function counts motifs of different classes
#It only counts motifs where the direction makes
#a difference, so where the main effects are different
#signs. It records whether the target or the source
#has the negative or posisive effect


count.motifs.directed <- function(motif.obj, plot.table = FALSE, write.table = FALSE, table.file = "Motifs.txt"){


	num.pheno <- length(motif.obj[[1]])
	motif.table <- matrix(NA, nrow = num.pheno, ncol = 4)
	
	colnames(motif.table) <- c("enhance-pos", "enhance-neg", "suppress-pos", "suppress-neg")
	rownames(motif.table) <- names(motif.obj[[1]])
	
	num.motifs <- rep(NA, num.pheno)
	names(num.motifs) <- names(motif.obj[[1]])
	for(ph in 1:num.pheno){
		
		main.diff <- which(motif.obj[[2]][[ph]][,2] != motif.obj[[2]][[ph]][,3])
		
		pos.int <- which(motif.obj[[2]][[ph]][,1] == 1)
		neg.int <- which(motif.obj[[2]][[ph]][,1] == -1)
		
		target.pos <- which(motif.obj[[2]][[ph]][,3] == 1)
		target.neg <- which(motif.obj[[2]][[ph]][,3] == -1)
		
		enhance.pos <- intersect(pos.int, target.pos)
		enhance.neg <- intersect(pos.int, target.neg)
		suppress.pos <- intersect(neg.int, target.pos)
		suppress.neg <- intersect(neg.int, target.neg)
		
		motif.table[ph, 1] <- length(intersect(main.diff, enhance.pos))
		motif.table[ph, 2] <- length(intersect(main.diff, enhance.neg))
		motif.table[ph, 3] <- length(intersect(main.diff, suppress.pos))
		motif.table[ph, 4] <- length(intersect(main.diff, suppress.neg))		

		num.motifs[ph] <- length(main.diff)		

		}	
	
	motif.prop <- motif.table/num.motifs
		
	if(plot.table){
		par(mfrow = c(1,2))
		imageWithText(motif.table, col.names = colnames(motif.table), row.names = rownames(motif.table), cex = 1)
		imageWithText(motif.prop, col.names = colnames(motif.table), row.names = rownames(motif.table), cex = 1)
		}
	
	
	result <- list(motif.table, motif.prop, num.motifs)
	names(result) <- c("motif.counts", "motif.proportions", "total.motifs")
	return(result)
	
	
	
}