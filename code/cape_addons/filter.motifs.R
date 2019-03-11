filter.motifs <- function(path = "."){
	
	motif.files <- get.files(want = c("Motifs", ".txt"))
	for(i in 1:length(motif.files)){
		motifs <- read.table(motif.files[i], sep = "\t", stringsAsFactors = FALSE, header = TRUE)
		
		main.same <- apply(motifs, 1, function(x) sign(as.numeric(x[4])) == sign(as.numeric(x[5])))
		if(length(which(main.same)) > 0){
			main.same.subtable <- motifs[which(main.same),,drop=FALSE]
			int.diff <- apply(main.same.subtable, 1, function(x) sign(as.numeric(x[4])) != sign(as.numeric(x[7])))
			if(length(int.diff) > 0){
				int.diff.subtable <- main.same.subtable[which(int.diff),,drop=FALSE]
				write.table(int.diff.subtable, paste0("main.int.diff.", motif.files[i]), sep = "\t", quote = FALSE, row.names = FALSE)
				}
			}
	}
	
	
	
}