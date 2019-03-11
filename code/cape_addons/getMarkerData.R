#getMarkerData.R
#Justin Hendrick, JustinJHendrick@gmail.com, Justin.Hendrick@jax.org
#8-9-13

#given a marker, return a named list with:
#    name, chromosome, cM, covar
getMarkerData = function(data.obj, marker, covar, covs, bp, bpTocM) {
    out = list(marker)
    names(out)[length(out)] = "name"
    
    #phenotypes or eigentraits?
    if(data.obj$transform.to.phenospace) {
        phenoNames = colnames(data.obj$pheno)
    } else {
        phenoNames = colnames(data.obj$ET)
    }
    
    markerI = which(data.obj$marker.names == marker)#find marker index
    chromo = as.numeric(data.obj$chromosome[markerI])
    if (length(chromo) == 0) {
        chromo = NULL
    } else if (chromo == 0) {#covar
        chromo = max(as.numeric(data.obj$chromosome)) + which(covs == marker)[1]
    }
    loc = data.obj$marker.location[markerI]
    
    if(length(loc) == 0) {
        cM = NULL
    } else {
        if(bp) { #convert bp to cM
            cM = loc * bpTocM
        } else {
            cM = loc
        }
        if(cM > 100) {
            cM = 100
        }
    }
    
    if(covar) {
        cM = 50
    }
    #only add non-null elements
    if(!(is.null(chromo) || is.na(chromo))) {
        out = append(out, chromo - 1)#silly biologists count from 1
        names(out)[length(out)] = "chromosome"
    }
    if(!is.null(cM)) {
        out = append(out, cM)
        names(out)[length(out)] = "cM"
    }
    if(covar) {
        out = append(out, which(covs == marker)[1])
        names(out)[length(out)] = "cov"
    }
    
    return(out)
}