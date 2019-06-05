#' Generate the data for simulation 3
#'
#' @param U number, effect size
#' @param seed integer, random seed
#'
#' @return a data list with following items
#'   Y: observation vector of length n;
#'   X: n by p covariate matrix;
#'   U1: p by m1 matrix, indicate feature membership
#'   U2: m1 by m2 matrix, indicate level-1 group membership
#'   types: the type of each feature, different type can have different s^2
#'   true_beta: true coefficients
#' @export
geneSimu3 <- function(U=0.2, seed=1){

  set.seed(seed)
  N <- 200         # number of samples
  P_dup <- 306     # number of features (including "duplicated" features)
  P <- 300
  m2 <- 10         # number of level-2 groups
  m1_dup <- 102    # number of level-1 groups (including "duplicated" groups)
  
  kl_dup <- c(10, 11, 10, 11, rep(10, 6))   # number of level-1 groups in each level-2 group
  pk_dup <- rep(3, m1_dup)   # number of features in each level-1 group
  pl_dup <- c(30, 33, 30, 33, rep(30, 6))   # number of features in each level-2 group
  X_dup <- matrix(NA, N, P_dup)  # feature matrix
  Y <- rep(NA, N)            # response
  l_index <- rep(seq(1, m2), times=pl_dup)       # level-2 index for each feature
  l_index_k <- rep(seq(1, m2), times=kl_dup)       # level-2 index for each level-1 group
  k_index_dup <- rep(seq(1, m1_dup), times=pk_dup)    # level-1 index for each feature
  
  true_c_ind <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0) # indicate if level-2 group contribute to response
  true_g_ind_dup <- c(c(rep(1, 6), rep(0, 4)), c(rep(1, 7), rep(0, 4)), 
    c(rep(1, 6), rep(0, 4)), c(rep(1, 7), rep(0, 4)), c(rep(1, 6), rep(0, 4)), 
                rep(0, 50))   # indicate if level-1 group contribute

  # directly simulate true beta (without duplication)
  true_beta_ind <- c(rep(c(rep(1,18), rep(0, 12)), times=5), 
                   rep(0, 150)) # indicate if feature contribute
  level_ind <- c(rep(c(rep(3, 12), rep(2, 6), rep(1, 12)), times=5), 
               rep(1, 150))     # effect size level of each feature

  true_beta <- rep(0, P)    
  true_beta[level_ind==1] <- 0
  true_beta[level_ind==2] <- sample(c(-1, 1), sum(level_ind==2), replace=TRUE, 
    prob=c(0.5, 0.5))*runif(sum(level_ind==2), U, 2*U)
  true_beta[level_ind==3] <- sample(c(-1, 1), sum(level_ind==3), replace=TRUE, 
    prob=c(0.5, 0.5))*runif(sum(level_ind==3), 2*U, 3*U)

  for(n in 1:N){
    zl <- rnorm(m2, 0, sqrt(0.2))
    zl_rep_dup <- rep(zl, times=pl_dup)
    zk_dup <- rnorm(m1_dup, 0, sqrt(0.3))
    zk_rep_dup <- rep(zk_dup, times=pk_dup)
    e_dup <- rnorm(P_dup, 0, sqrt(0.5))
    X_dup[n, ] <- zl_rep_dup + zk_rep_dup + e_dup
  }

	# set the duplicated feature to be same, and equals to average
  over_group_index1 <- c(1, 11) 
  over_group_index2 <- c(22, 32)
  X_dup[, k_index_dup==over_group_index1[1]] <- X_dup[, k_index_dup==over_group_index1[2]] <- 
    (X_dup[, k_index_dup==over_group_index1[1]]
    +X_dup[, k_index_dup==over_group_index1[2]])/2
  X_dup[, k_index_dup==over_group_index2[1]] <- X_dup[, k_index_dup==over_group_index2[2]] <- 
    (X_dup[, k_index_dup==over_group_index2[1]]
    +X_dup[, k_index_dup==over_group_index2[2]])/2
  
  # true feature matrix   
  X <- X_dup[, -c(which(k_index_dup %in% c(over_group_index1[2], 
    over_group_index2[2])))]
  Y <- X %*% true_beta + rnorm(N)
  P <- ncol(X)
  m1 <- m1_dup - 2
  pk <- rep(3, m1)
  kl <- rep(10, m2)

  # create feature group (U1) matrix
  U1 <- matrix(0, P, m1)
  for(k in 1:m1){
  	if(k==1){
  		U1[1:pk[1], 1] <- 1
  	}else{
  		U1[(sum(pk[1:(k-1)])+1):(sum(pk[1:k])), k] <- 1
  	}
  }

  U2 <- matrix(0, m1, m2)
  for(l in 1:m2){
  	if(l==1){
  		U2[1:kl[1], 1] <- 1
  	}else{
  		U2[(sum(kl[1:(l-1)])+1):(sum(kl[1:l])), l] <- 1
  	}
  }
  U2[1, 2] <- 1  # level-1 group 1 belongs to both level-2 group 1 and 2
  U2[21, 4] <- 1 # level-1 group 21 belongs to both level-2 group 3 and 4

	return(list(Y=Y, X=X, U1=U1, U2=U2, true_beta=true_beta))
}
