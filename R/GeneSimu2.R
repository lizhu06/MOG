#' Generate the data for simulation 2
#'
#' @param seed numeric number, random seed
#'
#' @return a data list with following items
#' 	 Y: observation vector of length n;
#'    X: n by p covariate matrix;
#'    U1: p by m1 matrix, indicate feature membership
#'    types: the type of each feature, different type can have different s^2
#'    true_beta: true coefficients
#' @export
#'
geneSimu2 <- function(seed=1){
	set.seed(seed)
	# We will create wo features that are shared by two groups. Each of these
	#		two features will have two "duplicated" features generated first. Then,
	#		the true feature will be the average of the "duplicated" features
	N <- 125        # number of samples
	P_dup <- 202    # number of features (including duplicated features)        
	m1 <- 10        # number of level-1 groups 
	pk_dup <- c(20,21,20,21, rep(20,6)) # number of features in each group (including duplicated features)  
	X_dup <- matrix(NA, N, P_dup) # feature matrix (including duplicated features)        
	Y <- rep(NA, N)
	true_group_ind <- c(rep(1,5), rep(0,5)) # indicate if group contribute to outcome
	true_beta_nonzero <- rnorm(50, 0, sqrt(5)) # sample all non-zero beta's (including duplicated features)
	true_beta <- c(c(true_beta_nonzero[1:20]),
	             c(true_beta_nonzero[21:30], rep(0,10)),
	             c(true_beta_nonzero[31:40], rep(0,10)),
	             c(true_beta_nonzero[41:45], rep(0,15)),
	             c(true_beta_nonzero[46:50], rep(0,15)), rep(0,100))
	for(n in 1:N){
	  z <- rnorm(m1, 0, 1)
	  e <- rnorm(P_dup, 0, 1)
	  for(k in 1:m1){
	    if(k == 1){
	      X_dup[n, 1:pk_dup[1]] <- (z[k] + e[1:pk_dup[1]])/sqrt(1+1)
	    }else{
	      X_dup[n, ((sum(pk_dup[1:(k-1)])+1):(sum(pk_dup[1:k])))] <- 
	      	(z[k]+e[((sum(pk_dup[1:(k-1)])+1):(sum(pk_dup[1:k])))])/sqrt(1+1)
	    }
	  }
	}

	# set the duplicated feature to be same, and equals to average
	overlap_index1 <- c(1, 21)  # duplicated feature 1
	overlap_index2 <- c(42, 62) # duplicated feature 2
	X_dup[, overlap_index1] <- cbind(apply(X_dup[, overlap_index1], 1 ,mean),
	  apply(X_dup[,overlap_index1],1,mean)) 
	X_dup[,overlap_index2] <- cbind(apply(X_dup[,overlap_index2],1,mean),
	  apply(X_dup[,overlap_index2],1,mean))
	#Y <- X_dup %*% true_beta_dup + rnorm(N)  # subject to change

	# true feature matrix 
	X <- X_dup[,-c(overlap_index1[2],overlap_index2[2])]
	Y <- X %*% true_beta + rnorm(N)  # subject to change
	P <- length(true_beta) # number of true features
	pk <- rep(20, 10)
	
	# create feature group (U1) matrix
	U1 <- matrix(0, P, m1)
	for(k in 1:m1){
		if(k==1){
			U1[1:pk[1], 1] <- 1
		}else{
			U1[(sum(pk[1:(k-1)])+1):(sum(pk[1:k])), k] <- 1
		}
	}
	U1[1, 2] <- 1 # feature 1 belongs to group 1 and 2
	U1[41, 4] <- 1 # feature 41 belongs to group 3 and 4
	types <- rep(1, P)
	return(list(Y=Y, X=X, U1=U1, types=types, true_beta=true_beta))
}
