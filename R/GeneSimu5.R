#' Generate the data for simulation 5
#'
#' @param N integer, sample size
#' @param s2 numeric number, specifying s2 in SOG model
#' @param m1_signal integer, number of level-1 groups that are predictive of outcome 
#' @param m1_noise integer, number of level-1 groups that are irrelevant to outcome 
#' @param p_perk integer, number of features inside each group
#' @param a numeric number, a in beta distribution (prior of pi_k^0)
#' @param b numeric number, b in beta distribution (prior of pi_k^0)
#' @param seed integer, specifying the random seed
#'
#' @return a data list with following items
#'   Y: observation vector of length n;
#'   X: n by p covariate matrix;
#'   U1: p by m1 matrix, indicate feature membership
#'   types: the type of each feature, different type can have different s^2
#'   true_beta: true coefficients
#'   pi0_k: true pi_k^0
#' @export
#'
geneSimu5 <- function(N=300, s2=5, m1_signal, m1_noise, p_perk=20, a, b, seed=123){
	set.seed(seed)

	m1 <- m1_signal + m1_noise
	pk <- c(rep(p_perk, m1_signal/2), rep(p_perk/4, m1_signal/2), rep(p_perk/2, m1_noise))
	cumsum_pk <- cumsum(pk)
	P <- sum(pk)
	U1 <- matrix(0, P, m1)
	for(k in 1:m1){
		if(k==1){
			U1[1:cumsum_pk[1], k] <- 1
		}else{
			U1[(cumsum_pk[k-1]+1):cumsum_pk[k], k] <- 1
		}	
	}
	
	# denote if it's signal or noise
	U1_signal <- matrix(0, P, m1)
	U1_signal <- U1[, 1:m1_signal]

	# generate X
	k_index <- unlist(lapply(1:length(pk), function(x) rep(x, pk[x])))
	X <- matrix(NA, N, P)
	Y <- rep(NA, N)

	for(n in 1:N){
	  z <- rnorm(m1, 0, 1)
	  e <- rnorm(P, 0, 1)
	  X[n, ] <- (z[k_index]+e)/sqrt(2)
	}

	# generate pi1 and pi0_k (prob of zero)
	pi1 <- 1 - m1_signal/m1
	pi0_k <- c(rbeta(m1_signal, a, b), rep(0.5, m1_noise))
	gamma_0 <- unlist(lapply(1:m1, function(k) rbinom(pk[k], 1, 1 - pi0_k[k])))

	# true_beta
	true_beta <- rep(0, P)
	true_beta[gamma_0==1] <- rnorm(sum(gamma_0), 0, s2)
	true_beta[k_index > m1_signal] <- 0
	# get back to non-duplicated X and beta
	Y <- X %*% true_beta + rnorm(N)

	types <- rep(1, P)
	return(list(Y=Y, X=X, U1=U1, types=types, true_beta=true_beta, pi0_k=pi0_k))
}
