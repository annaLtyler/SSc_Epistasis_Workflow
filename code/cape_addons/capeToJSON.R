#capeToJSON.R
#Justin Hendrick
#JustinJHendrick@gmail.com
#Justin.Hendrick@jax.org
#8-9-13
#convert cape RData file to JSON for vis
#output has: links, nodes, nCovar, phenos, mins and maxes of phenos and covars, chromoLens, pVal

capeToJSON = function(data.obj, outDir = ".", outFileName = "capeOut.json") {

	cur.dir <- getwd() #get the current working directory so we can return here after writing the output
	setwd(outDir)
	
    #load necessary code
    require(RJSONIO)
    #library(cape)
    #source("../../writeVariantInfluences.R")
    # source("genEffectPlotData.R") #adapted from plot.effects.R
    # source("get.interaction.error.R") #an addon not included in package
    # source("getMarkerData.R") #retrieves neccessary data about a specific marker
    whichBlock = function(data.obj, mIndex) {#find index of block that contains mIndex
        i = 0 #i is a js index. Starts at 0
        for(b in data.obj$linkage.blocks.collapsed) {
            for(n in b) {
                if(mIndex %in% which(colnames(data.obj$geno) == n)) {
                    return(i)
                }
            }
            i = i + 1
        }
        return(-1) #phenotypes and covars aren't in blocks
    }

    pVal = data.obj$network.p.or.q
    if(is.null(pVal)) {
        stop("get.network() must be run before converting to JSON")
	    }
    
    if(data.obj$transform.to.phenospace) {
        phenoNames = colnames(data.obj$pheno)
    } else {
        phenoNames = colnames(data.obj$ET)
    }
    #load the variant influences
    varInf = writeVariantInfluences(data.obj, write.file = FALSE, p.or.q = pVal, mark.covar = TRUE)
    rownames(varInf) = NULL
    colnames(varInf) = NULL

    #source of these two functions: kBroman
    cat0 <- function(...) cat(..., sep="", file=outFileName) #open new file
    catAppend <- function(...) cat(..., sep="", file=outFileName, append=TRUE) #append to already open file

    nCovar = 0
    covs = c()
    markers = c()
    isCov = c()
    bp = FALSE
    if(max(data.obj$marker.location) > 150) { #are marker positions in base pairs?
        bp = TRUE
        bpTocM = 100 / (max(data.obj$marker.location) - min(data.obj$marker.location))
    }
    
    nodeNames = c(data.obj$marker.names, phenoNames, names(data.obj$linkage.blocks.collapsed))
    cat0("{\"links\":[")
    for(r in 1:nrow(varInf)) {#for each row. write the link object
        sourceStr = varInf[r, 1]
        targetStr = varInf[r, 2]
        mSource = c()
        mTarget = c()
        
        sCov = FALSE
        tCov = FALSE
        if (substr(sourceStr, nchar(sourceStr), nchar(sourceStr)) == "*") {#last char is * means covar
            sourceStr = substr(sourceStr, 1, nchar(sourceStr) - 1) #remove *
            sCov = TRUE
        }
        if (substr(targetStr, nchar(targetStr), nchar(targetStr)) == "*") {#last char is * means covar
            targetStr = substr(targetStr, 1, nchar(targetStr) - 1) #remove *
            tCov = TRUE
        }

        sourceI = which(nodeNames == sourceStr)
        targetI = which(nodeNames == targetStr)
        if(!sourceI %in% markers) { #if not in, add it. This will be the list of marker indices with data on them
            markers = append(markers, sourceI)
            isCov = append(isCov, sCov)
        }
        mSource[2] = which(markers == sourceI)
        
        if(!targetI %in% markers) { #same for target
            markers = append(markers, targetI)
            isCov = append(isCov, tCov)
        }
        mTarget[2] = which(markers == targetI)
        
        sBlock = whichBlock(data.obj, sourceI) #which block is it in?
        tBlock = whichBlock(data.obj, targetI)
        if(sBlock != -1) { #in a block. put into mSource
            mSource[1] = sBlock
        } else { #-1 means it's not in a block. copy marker index
            mSource[1] = mSource[2] + length(data.obj$linkage.blocks.collapsed) - 1
        }
        if(tBlock != -1) { #same for target
            mTarget[1] = tBlock
        } else {
            mTarget[1] = mTarget[2] + length(data.obj$linkage.blocks.collapsed) - 1
        }

        eff  = as.numeric(varInf[r, 3])
        sig  = as.numeric(varInf[r, 5])
        qval = as.numeric(varInf[r, 7]) #p.or.q
        if(sourceStr %in% data.obj$marker.names && targetStr %in% data.obj$marker.names) {#both are markers, not phenos
            effPlot = genEffectPlotData(data.obj, sourceStr, targetStr)
            out = list(mSource, mTarget, eff, sig, qval, effPlot)
            names(out) = list("source", "target", "eff", "sig", "qval", "ePlot")
        } else {
            out = list(mSource, mTarget, eff, sig, qval)
            names(out) = list("source", "target", "eff", "sig", "qval")
        }
        out[["source"]][2] = out[["source"]][2] + length(data.obj$linkage.blocks.collapsed) - 1 #R is 1 indexed and js is 0 indexed
        out[["target"]][2] = out[["target"]][2] + length(data.obj$linkage.blocks.collapsed) - 1
        catAppend(toJSON(out))
        if (r != nrow(varInf)) {#not last one. put a comma
            catAppend(",")
        } else {#last one, end of array
            catAppend("],")
        }
    }
    
    #make vector of covariate(s)
    for(i in 1:length(markers)) {
        if (isCov[i]) {
            if(!markers[i] %in% covs) {#new cov
                nCovar = nCovar + 1
                covs[nCovar] = nodeNames[markers[i]]
            }
        }
    }
    
    #nCovar should be equal to the sum of the flags
    #this case arises when a covariate has no significant influences
    #and therefore, is not included in writeVariantInfluences table
    nShouldBe = sum(data.obj$covar.flags)
    if(nCovar != nShouldBe) {
        #if it isn't, make it so.
        nCovar = nShouldBe
    }
    
    #nodes: blocks then markers
    catAppend("\"nodes\":[")
    blocks = data.obj$linkage.blocks.collapsed
    names(blocks) = NULL
    
    for(b in blocks) {
        for(i in 1:length(b)) {
            b[i] = which(markers == which(colnames(data.obj$geno) == b[i])) + length(data.obj$linkage.blocks.collapsed) - 1;
        }
        catAppend(toJSON(b), ",")
    }
    
    for(i in 1:length(markers)) {
        catAppend(toJSON(getMarkerData(data.obj, nodeNames[markers[i]], isCov[i], covs, bp, bpTocM)))
        if(markers[length(markers)] != markers[i]) {#not last one
            catAppend(",")
        } else { #last one
            catAppend("],") #end of nodes
        }
    }
    
    #build chromoLens based on number of markers per chromosome. only relative size matters
    chromoLens = c()
    chromo = as.numeric(data.obj$chromosome)
    i = 1
    for(ch in chromo) {
        if(ch == 0 && nCovar > 0) {#covar
            ch = max(chromo) + i
            i = i + 1  #add one more for the next covar
        }
        
        if(ch > 0 && (is.null(chromoLens[ch]) || is.na(chromoLens[ch]))) {#no value yet counts as zero, make it one
            chromoLens[ch] = 1
        } else if (ch > 0) {
            chromoLens[ch] = chromoLens[ch] + 1      #add marker to this chromosome
        }
    }
    
    chromoLens[is.na(chromoLens)] = .1 #turn NAs to .1
    
    #geno min and max
    gMin = min(data.obj$geno[, -(ncol(data.obj$geno) - nCovar + 1):-ncol(data.obj$geno)], na.rm = TRUE)
    gMax = max(data.obj$geno[, -(ncol(data.obj$geno) - nCovar + 1):-ncol(data.obj$geno)], na.rm = TRUE)

    #covar mins and maxs
    if(nCovar > 0) {
        cMins = c()
        cMaxs = c()
        for(co in 1:nCovar) {
            cMins[co] = min(data.obj$geno[, ncol(data.obj$geno) - nCovar + co], na.rm = TRUE)
            cMaxs[co] = max(data.obj$geno[, ncol(data.obj$geno) - nCovar + co], na.rm = TRUE)
        }
    }
    
    catAppend("\"nCovar\":", nCovar, ",")
    catAppend("\"phenos\":", toJSON(phenoNames), ",")
    catAppend("\"gMin\":", gMin, ",")
    catAppend("\"gMax\":", gMax, ",")
    if(nCovar > 0) {
        catAppend("\"cMins\":", toJSON(cMins), ",")
        catAppend("\"cMaxs\":", toJSON(cMaxs), ",")
    }
    catAppend("\"chromoLens\":", toJSON(chromoLens), ",")
    catAppend("\"pVal\":", pVal)
    catAppend("}")
    
    setwd(cur.dir)
}