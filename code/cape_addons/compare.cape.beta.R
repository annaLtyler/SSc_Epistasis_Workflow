#This function compares cape coefficients to standard epistasis
#coefficients. The function requires a cape data object and 
#a cross object generated using calc.p.beta using the same
#data and arguments


compare.cape.beta <- function(cape.obj, beta.obj){
	
	cape.influences <- cape.obj$var.to.var.influences
	beta.terms <- beta.obj$all.var.to.var.p.val
	
	std.m12 <- cape.influences[,"m12"]/cape.influences[,"m12.std.dev"]
	std.m21 <- cape.influences[,"m21"]/cape.influences[,"m21.std.dev"]
	
	for(i in 1:length(beta.terms)){
		pheno.betas <- beta.terms[[i]][,"Effect"]/beta.terms[[i]][,"SE"]
		quartz(width = 10, height = 5);par(mfrow = c(1,2))
		plot(std.m12, pheno.betas, main = names(beta.terms)[i])
		abline(v = 0); abline(h = 0)
		plot(std.m21, pheno.betas, main = names(beta.terms)[i])
		abline(v = 0); abline(h = 0)
		}
	
	
	
	
}