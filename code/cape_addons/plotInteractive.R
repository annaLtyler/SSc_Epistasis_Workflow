#plotInteractive.R
#Justin Hendrick, JustinJHendrick@gmail.com, justin.hendrick@jax.org
#8-9-13

#optionally writes capeOut.json, then starts a python server and opens it in the default browser
#plotInteractive.js and index.html should be in the working directory

plotInteractive = function(data.obj, outDir = ".", pythonDir = "", port = 8888, tutorial = FALSE) {
    require(RCurl)
    
    outFileName = "capeOut.json"
    
      if(length(unique(data.obj$geno[,1])) > 3){
    	message("I am detecting more than 3 genotypes per marker. Do you want to bin the genotypes to 0, 0.5 and 1 values (y/n)?\n")
    	bin <- readLines(n = 1)
    	if(bin == "y"){
    		cat("Binning Genotypes...\n")
    		data.obj <- bin.geno(data.obj)
    		}else{
    		cat("Not Binning Genotypes\n")
    		message("Not binning the genotypes may cause this script to take a long time and give uninterpretable effect plots.\n")
    		message("You can bin markers using bin.geno()\n")
    		}
    	}
 
    
    files = list.files(path = outDir)
    if(!"index.html" %in% files) {
        stop("index.html is missing")
    	}

    if(!"plotInteractive.js" %in% files) {
        stop("plotInteractive.js is missing")
    	}
    
      
    cat("writing json...\n")
    capeToJSON(data.obj, outDir = outDir, outFileName = outFileName)
	cat(" done\n")
    
    files = list.files(path = outDir)
    if(!outFileName %in% files) {
        stop(paste(outFileName, "is missing"))
 	   }
    
   
    if(!url.exists(paste0("http://localhost:", port))) { #only create a new server if the port is empty
        if(nchar(pythonDir) > 0 && substr(pythonDir, nchar(pythonDir), nchar(pythonDir)) != "/") { #change /path/foo/bar to /path/foo/bar/
            pythonDir = paste0(pythonDir, "/")
        	}
        
        #call python -V and store the result
        vString = system(paste0(pythonDir, "python -V"), intern = TRUE)
        #extract major version number from string such as this "Python 3.3.2"
        pythonV = as.numeric(strsplit(strsplit(vString, " ", TRUE)[[1]][2], ".", TRUE)[[1]][1])
        
        cat(paste("starting http server with python on port", port, "\n"))
        if(!is.na(pythonV) && pythonV >= 3) { #python 3 and up
            system(paste0(pythonDir, "python -m http.server ", port), wait = FALSE)
        } else if(!is.na(pythonV) && pythonV < 3) { #python 2 and down
            system(paste0(pythonDir, "python -m SimpleHTTPServer ", port), wait = FALSE)
        } else {
            stop(paste("Unknown python version", pythonV))
        }
        Sys.sleep(2) #wait for the server to start up because it was called asynchronously
    } else {
        cat(paste("server already exists on port", port, "Opening it\n"))
    }
    
    if(tutorial) {
        browseURL(paste0("http://localhost:", port, "/tutorial"))
    } else {
        browseURL(paste0("http://localhost:", port))
    }
}