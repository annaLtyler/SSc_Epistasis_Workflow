#This function reads

read.GTEx.output <- function(filename, tissues.of.interest = NULL){

	data.table <- read.csv(filename, stringsAsFactors = FALSE)
	
	locate.tissue <- function(tissue.name){
		tissue.idx <- grep(tissue.name, data.table[,"Tissue"], ignore.case = TRUE)
		tissue.table <- cbind(tissue.idx, rep(tissue.name, length(tissue.idx)))
		return(tissue.table)
		}
	
	if(!is.null(tissues.of.interest)){
		tissue.locale <- Reduce("rbind", lapply(tissues.of.interest, function(x) locate.tissue(x)))
		
		if(length(tissue.locale) > 0){
			sub.table <- data.table[as.numeric(tissue.locale[,1]),,drop=FALSE]
			sub.table[,"Tissue"] <- tissue.locale[,2]
			}else{
			return("No genes affected in the tissues of interest")	
			}
		}else{
		sub.table <- data.table	
		}
	
	u_genes <- unique(sub.table[,"Gene.Symbol"])
	tissue.list <- lapply(u_genes, function(x) sub.table[which(sub.table[,"Gene.Symbol"] == x),"Tissue"])
	tissue.text <- sapply(tissue.list, function(x) paste(unique(x), collapse = ";"))
	final.table <- cbind(rep(data.table[1,"SNP.Id"],length(u_genes)), u_genes, tissue.text)
	colnames(final.table) <- c("SNP", "gene", "tissue")
	return(final.table)
	}