#This script selects phenotypes that have correlation coeffieicnts
#within a given range

select.correlations <- function(data.obj, min.cor = 0.4, max.cor = 0.8, label = NULL){
	

	cor.table <- cor(data.obj$pheno, use = "complete")

	if(is.null(min.cor)){
		sel.cor <- which(abs(cor.table) <= max.cor, arr.ind = TRUE)
		pdf.label <- paste("Greater.Than", min.cor, ".", label, sep = "")
		}

	if(is.null(max.cor)){
		sel.cor <- which(abs(cor.table) >= min.cor, arr.ind = TRUE)
		pdf.label <- paste("Less.Than.", max.cor, ".", label, sep = "")
		}
		
	if(!is.null(min.cor) && !is.null(max.cor)){
		sel.cor <- which(abs(cor.table) >= min.cor & abs(cor.table) <= max.cor, arr.ind = TRUE)
		pdf.label <- paste("Between.", min.cor, ".and.", max.cor, ".", label, sep = "")
		}

	
	
	get.cor <- function(pheno.pair){
		pheno1 <- rownames(cor.table)[pheno.pair[1]]
		pheno2 <- rownames(cor.table)[pheno.pair[2]]
		return(c(pheno1, pheno2, cor.table[pheno.pair[1], pheno.pair[2]]))
		}
	
	final.table <- t(apply(sel.cor, 1, get.cor))
	
	write.table(final.table, paste("Correlations.", pdf.label, ".txt", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")	
	
	included.phenotypes <- unique(c(final.table[,1], final.table[,2]))
	invisible(included.phenotypes)
	
	}