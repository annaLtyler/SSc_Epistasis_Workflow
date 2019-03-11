#This function plots genes in position on a chromosome.
#gene.table must have the following columns: external_gene_name, 
#chromosome_name, and start_position or end_position

plot.gene.pos <- function(gene.table, highlight.genes = NULL, highlight.col = "red", plot.label = NULL){
	
	u_chr <- unique(gene.table[,"chromosome_name"])
	sub.tables <- lapply(u_chr, function(x) gene.table[which(gene.table[,"chromosome_name"] == x),])
	
	
	for(i in 1:length(sub.tables)){
		pos.col <- grep("position", colnames(sub.tables[[i]]))[1] #take either start or stop
		chr.min <- min(sub.tables[[i]][,pos.col])
		chr.max <- max(sub.tables[[i]][,pos.col])
		colV <- rep("gray", nrow(sub.tables[[i]]))
			
		if(!is.null(highlight.genes)){
			highlight.pos <- match(highlight.genes, sub.tables[[i]][,"external_gene_name"])
			colV[highlight.pos] <- highlight.col
			}
		
		plot.new()
		plot.window(xlim = c(chr.min, chr.max), ylim = c(0,1))
		points(x = sub.tables[[i]][,pos.col], y = runif(nrow(sub.tables[[i]]), 0,1), col = colV, pch = 16)
		axis(1)		
		mtext(paste(plot.label, "\nChr", u_chr[i]))
		
		}
	
	
}