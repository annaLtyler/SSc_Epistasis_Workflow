	get.delta <- function(markers,beta.m,se,beta.cov) {

	require("qpcR")
	require("corpcor")

		beta.main <- beta.m[,1:2]
		beta.inter <- beta.m[,3]
		n.rows <- dim(beta.main)[1] #No of ETs
		n.cols <- dim(beta.main)[2]

		se.main <- se[,1:2]
		se.inter <- se[,3]

		
		if(n.rows == n.cols){
			act_delta <- solve(beta.main)%*%beta.inter
			}else{
			act_delta <- pseudoinverse(beta.main)%*%beta.inter
			}

		return(act_delta)

}