#This function checks for markers that can't be used
#and directs the user to delete them using remove.unused.markers()

check.bad.markers <- function(data.obj, geno.obj = NULL){
	
	
	gene <- get.geno(data.obj, geno.obj)
	
	#we do not scan markers on the sex chromosomes
	#take these out here.
	x.locale <- grep("X", data.obj$chromosome, ignore.case = TRUE)
	y.locale <- grep("Y", data.obj$chromosome, ignore.case = TRUE)

	if(length(x.locale) > 0 || length(y.locale) > 0){
		stop("Sex chromosomes are not scanned in cape.\nPlease use remove.unused.markers() to remove these before proceeding.")
		}
		
	
	#take out markers with missing data
	# num.allele <- apply(gene, locus.dim, function(x) colSums(x, na.rm = TRUE))
	num.allele <- apply(gene, 3, function(x) sum(x, na.rm = TRUE))
	mono.allele <- which(num.allele == 0)
	if(length(mono.allele) > 0){
		stop("Some markers are invariant.\nPlease use remove.unused.markers() to remove these before proceeding.")
		}

	
	
	
	
	
}