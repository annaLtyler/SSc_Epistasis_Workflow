#This script adds a "change in" phenotype
#for any base phenotype measured twice.

add.change.in <- function(cross.obj, base.pheno){
	
	orig.pheno <- cross.obj$pheno
	orig.names <- colnames(cross.obj$pheno)
	
	#make sure the change in phenotype is not already present.
	change.name <- paste("change_in_", base.pheno, sep = "")
	already.present <- get.col.num(orig.pheno, change.name)
	
	if(length(already.present) > 0){
		message("The phenotype ", change.name, " already exists")
		return(cross.obj)
		}else{ #if the phenotype isn't here already, we must calculate it and add it.
			cat("Adding ", change.name, "\n")
			#first check to see that there were two measurements of the base pheno
			measurements <- grep(base.pheno, orig.names, ignore.case = TRUE)
			if(length(measurements) != 2){
				message("This base name does not have two measurements.")
				cat(orig.names[measurements],"\n")
				return(cross.obj)
				}else{ #if we do have two measurements, we can take the difference.
					measure1.name <- orig.names[measurements][grep(1, orig.names[measurements])]
					measure2.name <- orig.names[measurements][grep(2, orig.names[measurements])]
					measure1.locale <- grep(measure1.name, orig.names)
					measure2.locale <- grep(measure2.name, orig.names)
					#The difference is measurement 1 minus measurement 2
					new.pheno <- as.numeric(orig.pheno[,measure1.locale]) - as.numeric(orig.pheno[,measure2.locale])
					new.pheno.mat <- cbind(orig.pheno, new.pheno)
					colnames(new.pheno.mat) <- c(orig.names, change.name)
					cross.obj$pheno <- new.pheno.mat
					return(cross.obj)
					} #end case for having two measurements of the base phenotype
			
			} #end case for when the change phenotype is not already present
	

	} #end function
