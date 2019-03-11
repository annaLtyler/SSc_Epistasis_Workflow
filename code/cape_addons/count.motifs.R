#This function counts motifs of different classes


count.motifs <- function(motif.obj, plot.table = FALSE, write.table = FALSE, table.file = "Motifs.txt"){


	num.pheno <- length(motif.obj[[1]])
	motif.table <- motif.prop <- matrix(NA, nrow = num.pheno, ncol = 10)
	colnames(motif.table) <- colnames(motif.prop) <- c("enhancing", "enhancing-main-same", "enhancing-main-diff", "enhancing-main-pos", "enhancing-main-neg", "suppressive", "suppressive-main-same", "suppressive-main-diff", "suppressive-main-pos", "suppressive-main-neg")
	rownames(motif.table) <- rownames(motif.prop) <- names(motif.obj[[1]])
		
	final.table.to.write <- NULL
	num.motifs <- rep(NA, num.pheno)
	names(num.motifs) <- names(motif.obj[[1]])
	for(ph in 1:num.pheno){

		total.motifs <- length(motif.obj[[2]][[ph]][,1])
		
		table.to.write <- matrix(NA, nrow = total.motifs, ncol = 5)
		colnames(table.to.write) <- c("source", "target", "pheno", "interaction", "main effects")
		
		enhancing <- which(motif.obj[[2]][[ph]][,1] == 1)
		suppressing <- which(motif.obj[[2]][[ph]][,1] == -1)
		main.same <- which(motif.obj[[2]][[ph]][,2] == motif.obj[[2]][[ph]][,3])
		main.diff <- which(motif.obj[[2]][[ph]][,2] != motif.obj[[2]][[ph]][,3])
		main.pos <- intersect(which(motif.obj[[2]][[ph]][,2] == 1), which(motif.obj[[2]][[ph]][,3] == 1))
		main.neg <- intersect(which(motif.obj[[2]][[ph]][,2] == -1), which(motif.obj[[2]][[ph]][,3] == -1))
		en.main.same <- intersect(enhancing, main.same)
		supp.main.same <- intersect(suppressing, main.same)
		en.main.pos <- intersect(enhancing, main.pos)
		en.main.neg <- intersect(enhancing, main.neg)
		supp.main.pos <- intersect(suppressing, main.pos)
		supp.main.neg <- intersect(suppressing, main.neg)
		en.main.diff <- intersect(enhancing, main.diff)
		supp.main.diff <- intersect(suppressing, main.diff)

		if(length(en.main.pos) > 0){
			table.to.write[en.main.pos,1:3] <- motif.obj[[1]][[ph]][en.main.pos,]
			table.to.write[en.main.pos,4] <- "enhancing"
			table.to.write[en.main.pos,5] <- "positive"
			}

		if(length(en.main.neg) > 0){
			table.to.write[en.main.neg,1:3] <- motif.obj[[1]][[ph]][en.main.neg,]
			table.to.write[en.main.neg,4] <- "enhancing"
			table.to.write[en.main.neg,5] <- "negative"
			}

		if(length(en.main.diff) > 0){
			table.to.write[en.main.diff,1:3] <- motif.obj[[1]][[ph]][en.main.diff,]
			table.to.write[en.main.diff,4] <- "enhancing"
			table.to.write[en.main.diff,5] <- "different"
			}


		if(length(supp.main.pos) > 0){
			table.to.write[supp.main.pos,1:3] <- motif.obj[[1]][[ph]][supp.main.pos,]
			table.to.write[supp.main.pos,4] <- "suppressing"
			table.to.write[supp.main.pos,5] <- "positive"
			}

		if(length(supp.main.neg) > 0){
			table.to.write[supp.main.neg,1:3] <- motif.obj[[1]][[ph]][supp.main.neg,]
			table.to.write[supp.main.neg,4] <- "suppressing"
			table.to.write[supp.main.neg,5] <- "negative"
			}

		if(length(supp.main.diff) > 0){
			table.to.write[supp.main.diff,1:3] <- motif.obj[[1]][[ph]][supp.main.diff,]
			table.to.write[supp.main.diff,4] <- "suppressing"
			table.to.write[supp.main.diff,5] <- "different"
			}
		
		final.table.to.write <- rbind(final.table.to.write, table.to.write)

		motif.table[ph, 1] <- length(enhancing)
		motif.table[ph, 2] <- length(en.main.same)
		motif.table[ph, 3] <- length(en.main.diff)
		motif.table[ph, 4] <- length(en.main.pos)
		motif.table[ph, 5] <- length(en.main.neg)
		motif.table[ph, 6] <- length(suppressing)
		motif.table[ph, 7] <- length(supp.main.same)
		motif.table[ph, 8] <- length(supp.main.diff)
		motif.table[ph, 9] <- length(supp.main.pos)
		motif.table[ph, 10] <- length(supp.main.neg)
		
		motif.prop[ph,] <- motif.table[ph,]/total.motifs
		num.motifs[ph] <- total.motifs
		}	
	
	
	if(plot.table){
		imageWithText(motif.table, col.names = colnames(motif.table), row.names = rownames(motif.table), cex = 1)
		}
	if(write.table){
		write.table(final.table.to.write, table.file, quote = FALSE, sep = "\t", row.names = FALSE)
		}
		
	result <- list(motif.table, motif.prop, num.motifs)
	names(result) <- c("motif.counts", "motif.proportions", "total.motifs")
	return(list(result, final.table.to.write))
	
	
	
}