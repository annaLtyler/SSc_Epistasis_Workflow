#This function creates a multiple cross object
#given the path names, and the names of the objects
#to concatenate, it builds one object with marker
#information across all crosses, as well as the
#individual results specified~/Documents/Data/Frankel/Combined_Cross/Results/Just.RQ.B.to.H
# paths <- c("~/Documents/Data/Frankel/Combined_Cross/Results/Combined.Cross.no.duplicate.days.norm.B.to.H","~/Documents/Data/Frankel/Combined_Cross/Results/Just.Gria4","~/Documents/Data/Frankel/Combined_Cross/Results/Just.RQ.B.to.H", "~/Documents/Data/Frankel/Combined_Cross/Results/Just.Scn8a")
# internal.obj <- "max.var.to.pheno.influence"


concatenate.crosses <- function(paths, internal.obj){
	
	cross.names <- sapply(strsplit(paths, "/"), function(x) x[length(x)])
	concat.cross <- vector("list", length(paths)); names(concat.cross) <- cross.names
	all.marker.info <- vector("list", length(paths)); names(all.marker.info) <- cross.names

	base.dir <- getwd()
	for(i in 1:length(paths)){
		cat(cross.names[i], "\n")
		setwd(paths[i])
		cross <- readRDS("cross.RData")
		marker.info <- list(cross$marker.names, cross$chromosome, cross$marker.location)
		names(marker.info) <- c("marker.names", "chromosome", "marker.location")
		all.marker.info[[i]] <- marker.info
		sub.obj.locale <- which(names(cross) == internal.obj)
		if(length(sub.obj.locale) == 0){
			stop("I can't find the object in this cross. Please check the spelling.")
			}
		sub.obj <- cross[[sub.obj.locale]]
		concat.cross[[i]] <- sub.obj
		}
	
	#make sure we get all the markers from all crosses
	all.markers <- NULL
	for(i in 1:length(cross.names)){
		all.markers <- c(all.markers, all.marker.info[[i]]$marker.name)
		}
	all.markers <- unique(all.markers)
	
	all.chr <- NULL
	all.marker.location <- NULL
	#go through all markers and get chromosome and location information about them
	for(i in 1:length(all.markers)){
		#figure out which crosses contain the marker
		marker.in.cross <- sapply(all.marker.info, function(x) which(x$marker.names == all.markers[i])) 
		#pick one of them to use for the marker information
		marker.locale <- min(which(length(marker.in.cross) > 0))
		all.chr <- c(all.chr, all.marker.info[[marker.locale]]$chromosome[as.numeric(marker.in.cross[marker.locale])])
		all.marker.location <- c(all.marker.location, all.marker.info[[marker.locale]]$marker.location[as.numeric(marker.in.cross[marker.locale])])
		}

	concat.cross[["marker.names"]] <- all.markers
	concat.cross[["chromosome"]] <- all.chr
	concat.cross[["marker.location"]] <- all.marker.location	
	
	setwd(base.dir)
		
	return(concat.cross)
}