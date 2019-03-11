#Justin Hendrick, JustinJHendrick@gmail.com, Justin.Hendrick@jax.org, 8-9-13
#
#This function takes in the data object and two marker names
#and outputs the data necessary to draw effect plots
#
#this function started from plot.effects.R
#I removed stuff that I don't need then returning what would have been plotted

genEffectPlotData <- function(data.obj, marker, marker2) {
	error.type = "se"
	
    markers = vector()
	markers[1] <- colnames(data.obj$geno)[which(data.obj$marker.names == marker)]
    markers[2] <- colnames(data.obj$geno)[which(data.obj$marker.names == marker2)]
	marker.names <- c(marker, marker2)
	
    all.pheno.mat <- data.obj$pheno
	all.pheno <- colnames(all.pheno.mat)
    
    cells = list()
    for(ph in 1:length(all.pheno)) {
        errors <- get.interaction.error(data.obj$geno[,markers[1]], data.obj$geno[,markers[2]], all.pheno.mat[,ph], error.type = error.type)
        ylim <- c(min((errors$means - errors$se), na.rm = TRUE), max((errors$means + errors$se), na.rm = TRUE))
        #dat = tapply(all.pheno.mat[,ph], list(data.obj$geno[,markers[1]], data.obj$geno[,markers[2]]), mean)
        #rownames(dat) = NULL
        #colnames(dat) = NULL

        rownames(errors$means) = NULL
        colnames(errors$means) = NULL
        rownames(errors$se) = NULL
        colnames(errors$se) = NULL #removing names makes toJSON output arrays instead of objects.
        cells[[ph]] = list(errors$means, errors$se, ylim)
        names(cells[[ph]]) = list("dat", "errors", "ylim")
    }
    
    #find number of points that are not NA. Minimum of all three data sources
    n = min(sum(!is.na(data.obj$geno[, markers[1]])), sum(!is.na(data.obj$geno[, markers[2]])), sum(!is.na(all.pheno.mat[,ph])))
    out = list(n, cells)
    names(out) = list("n", "cells")
    return(out)
}