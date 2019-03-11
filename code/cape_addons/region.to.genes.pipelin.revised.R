#This is a pipeline that goes from a single Source-Target
#relationship to candidate genes

#========================================================================================
#load the data object you want to examine
#========================================================================================

region.to.genes.pipeline.revised <- function(data.obj, num.cand.interactions = 30, r2.thresh = 0.8){
	
	top.interactor.info(data.obj, num.pairs = num.cand.interactions, r2.thresh = r2.thresh)
	cat(paste("There is new output: Top.", num.cand.interactions, ".Interactions.r2.thresh.", r2.thresh, ".txt", sep = ""), "\n")	
	
	
	cat("\nTo translate marker regions to bp we need the bp location of each marker.\n")
	cat("This file should be named marker_bp.txt.\n")
	all.set <- readline(prompt = "Press enter when ready (x to quit):\n")
	
	if(all.set == "x"){
		return()
		}
	
	write.block.coord(data.obj, r2.thresh = r2.thresh, output.file = "Block.Coord.txt")

	cat("Now you can enter the block regions into MouseMine\n")
	cat("Look for genes in these regions with resonable annotations, and put candidates into IMP.")
	
}