#This function makes a table of linkage blocks with significant 
#main effects, it indicates the chromosome, the start and stop 
#position, the maximum effect size in the block, and the number 
#of genes in the interval.


main.effect.table <- function(data.obj, singlescan.obj){
	
	
	blocks <- data.obj$linkage.blocks.collapsed
	sig.levels <- unlist(singlescan.obj$alpha.thresh)
	
	singlescan.results <- singlescan.obj$singlescan.results
	
	
	get.sig.level <- function(max.effect){
		sig.locale <- which(sig.levels <= abs(max.effect))
		if(length(sig.locale) == 0){
			return("not significant")
			}else{
			return(names(sig.levels)[min(sig.locale)])	
			}	
		}
	
	get.block.stats <- function(block, trait.num){
		is.char <- as.logical(is.na(suppressWarnings(as.numeric(block[1]))))
		if(is.char){
			block <- get.marker.num(data.obj, block)
			}
		chr <- get.marker.chr(data.obj, block)
		pos <- get.marker.location(data.obj, block)
		start.pos <- min(pos)
		end.pos <- max(pos)
		ind.results <- singlescan.results[[trait.num]]
		marker.locale <- which(rownames(ind.results) %in% block)
		max.effect.locale <- which(ind.results[marker.locale,"t.stat"] == max(ind.results[marker.locale,"t.stat"]))
		max.effect <- ind.results[marker.locale[max.effect.locale],"slope"]/ind.results[marker.locale[max.effect.locale],"se"]
		max.marker <- rownames(ind.results)[marker.locale[max.effect.locale]]
		max.marker.name <- get.marker.name(data.obj, max.marker)
		sig.level <- get.sig.level(max.effect)
		table.row <- c(names(singlescan.results)[trait.num], length(block), unique(chr), start.pos, end.pos, max.marker.name, max.effect, sig.level)
		return(table.row)
		}	
		
	table.list <- vector(mode = "list", length = length(singlescan.results))
	for(i in 1:length(singlescan.results)){
		registerDoParallel()
		trait.table <- foreach(b = blocks, .combine = "rbind") %dopar% {
			get.block.stats(b, i)
			}
		colnames(trait.table) <- c("trait", "num.markers", "chr", "start.pos", "end.pos", "max.marker", "max.effect", "sig.level")
		rownames(trait.table) <- names(blocks)
		table.list[[i]] <- trait.table
		}
	
	names(table.list) <- names(singlescan.results)
	return(table.list)
	
}

