#This function reads in an OBO file
#and creates a list that I will parse with
#my own functions, since ontoCAT isn't cutting
#it

read.obo <- function(file.name){
	
	obo.table <- read.table(file.name, sep = "\t", stringsAsFactors = FALSE)
	term.locale <- which(obo.table[,1] == "[Term]")
	term.list <- vector(mode = "list", length = length(term.locale))

	last.term <- length(term.locale)

	for(i in 1:length(term.locale)){
		report.progress(i, length(term.locale))
		term.start <- term.locale[i]+1
		if(i < last.term){
			term.end <- term.locale[(i+1)] - 1
			}else{
			term.end <- nrow(obo.table)	
			}
		ind.term <- obo.table[term.start:term.end,1,drop=FALSE]
		split.terms <- strsplit(ind.term[,1], " ")
		att.names <- unlist(lapply(split.terms, function(x) x[1]))
		att.vals <- unlist(lapply(split.terms, function(x) paste(x[2:length(x)], collapse = " ")))
		att.table <- cbind(att.names, att.vals)
		term.list[[i]] <- att.table		
#		term.list[[i]] <- att.vals[1]
#		for(j in 2:length(att.names)){
#			attr(term.list[[i]], att.names[j]) <- att.vals[j]
#			}
		}	
	
	saveRDS(term.list, paste0(gsub("txt", "rds", basename(file.name))))	
	invisible(term.list)
}