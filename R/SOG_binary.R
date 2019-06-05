#' SOG function for binary outcome
#'
#' This function applies SOG to binary outcome
#'
#' @param X n by p feature matrix, where n is the number of samples, 
#' p is number of features
#' @param Y outcome vector of length n
#' @param U1 p by m1 matirx, denoting group membership of each feature
#' @param types numeric vector of length p, type index of each feature, allows different variance 
#' in priors of coefficient of different types 
#' @param lassoInit TRUE/FALSE, if TRUE, use lasso estimates as MCMC initials
#' @param pi2_weighted TRUE/FALSE, if TRUE, the prior of gamma_jk^0 will be weighted by R_j
#' @param pi2_prop_n numeric number, specifying n_0 in the proposal distribution of pi_k^0
#' @param MH_ind 1/0, if 1 MH is used to update pi_k^0; otherwise, Gibbs sampling is used
#' @param burnInIter number of iterations as burn-in period
#' @param keepIter number of iterations to be stored
#' @param maxIter maximum number of iterations
#' @param printInt print progress message every print_int iterations
#' @param seed integer, specifying the random seed
#' @param alpha2 numeric number, specifying alpha in the BETA distribution (prior of pi_k^0)
#' @param beta2 numeric number, specifying beta in the BETA distribution (prior of pi_k^0)
#' @return A list of MCMC output, together with feature selection probability 
#' and feature Bayesian qvalue
#' @export
SOG_binary <- function(X, Y, U1, types=NULL, 
  lassoInit=FALSE, pi2_weighted=FALSE, pi2_prop_n=10, MH_ind=1,
  burnInIter=1000, keepIter=3000, maxIter=10000, 
  printInt=1000,  seed=123,
  alpha2=1, beta2=1){
  
  set.seed(seed)
  s2Weighted <- FALSE
  
  # summarize the data
  P <- ncol(X)
  N <- length(Y)
  if(is.null(types)){
    types <- rep(1, P)
  }
  m1 <- ncol(U1)
  pk <- colSums(U1)   # number of features belonging to each group
  k_index <- unlist(sapply(1:P, function(j) which(U1[j,]==1)))

  U1_rowsum <- rowSums(U1)
  if(sum(U1_rowsum==0)>1){
    stop("Feature must belongs to at least one group, make a singleton group for feature that does not belong sto any groups")
  }
  overlap_index <- !(all(U1_rowsum==1))
  if(overlap_index){
    cat("Data include overlapping groups. \n")
  }
  if(overlap_index && pi2_weighted){
    if(MH_ind==0){
      stop("When overlapping groups exist and Bernoulli needs to be corrected, MH_ind have to be set 1. \n")
    }
  }
  # duplicate features if group overlaps
  w_weight <- rep(1, P)
  if(overlap_index){
    X_raw <- X
    feature_dup_index <- unlist(lapply(1:nrow(U1), function(j) rep(i, U1_rowsum[j])))
    k_index <- unlist(lapply(1:nrow(U1), function(i) which(U1[i,]==1)))
    X <- X_raw[, feature_dup_index]
    P <- ncol(X)
    w_weight <- rep(1, P)
    if(s2Weighted){
      w_weight <- unlist(lapply(1:length(U1_rowsum), function(j) rep(1/U1_rowsum[j], U1_rowsum[j])))
    }
    types <- types[feature_dup_index]
  }

  # weight for Bernoulli prior
  if(pi2_weighted){
    v_weight <- unlist(lapply(1:length(U1_rowsum), 
      function(j) rep(1/U1_rowsum[j], U1_rowsum[j])))
  }else{
    v_weight <- rep(1, ncol(X))
  }

  ##### hyper-parameters #####
  alpha1 <- 1
  beta1 <- 1
  #alpha2 <- 1
  #beta2 <- 1
  alpha_s <- 0
  beta_s <- 0
  alpha_sigma <- 0
  beta_sigma <- 0

  X <- cbind(rep(1, nrow(X)), X) # add intercept
  k_index <- c(m1+1, k_index) # intercept have fake group index
  types <- c(max(types)+1, types) # intercept has fake type

  ## Lasso 
  if(lassoInit == TRUE){
    fitcv <- glmnet::cv.glmnet(X[,-1], Y, family="binomial", 
      alpha=1, intercept=TRUE)
    coef <- coef(fitcv, s = "lambda.min")
    b <- coef[,1] # including intercept
    gamma_1 <- sapply(1:m1, function(x) max((b!=0)[k_index == x]))
    gamma_1_rep_per_feature <- c(1, gamma_1[k_index[-1]]) # including intercept
    gamma_2 <- b!=0   # include intercept
    pi1 <- max(sum(gamma_1)/length(gamma_1),0.1)
    pi2 <- sapply(1:m1, function(k) max(sum(gamma_2[k_index == k])/
      length(gamma_2[k_index == k]),0.1))
    pi2_rep_per_feature <- c(1, pi2[k_index[-1]])

    }else{
      b <- rnorm(P+1)
      pi1 <- 0.5
      pi2 <- rep(0.5, m1)
      pi2_rep_per_feature <- c(1,pi2[k_index[-1]])
      gamma_1 <- rbinom(m1,1,pi1)
      gamma_1_rep_per_feature <- c(1, gamma_1[k_index[-1]])
      gamma_2 <- rbinom(P+1,1,pi2_rep_per_feature)
  }
    
  sigma2 <- 1
  uni_types <- seq(1, max(types))
  s2 <- c(rep(1,length(uni_types)-1),100)  # variance vector for coefficients from different types

  beta <- gamma_1_rep_per_feature*gamma_2*b

  #s2Weighted_int <- s2Weighted*1
  pi2_weighted_int <- pi2_weighted*1
  #PreEstBeta_int <- PreEstBeta * 1

  get_pi2_logl <- function(k){
    sel_index <- which(k_index==k)
    return(sum(gamma_2[sel_index]*log(v_weight[sel_index]*pi2[k])) + 
      sum((1-gamma_2[sel_index])*log(1-v_weight[sel_index]*pi2[k])))
  }
  pi2_loglikeli <- sapply(1:m1, function(k) get_pi2_logl(k))

  res <- MCMC_sb(seed, burnInIter, keepIter, 
    m1, P, printInt, N, 
    alpha1, alpha2, alpha_sigma, alpha_s, 
    beta1, beta2, beta_sigma, beta_s, pi1, sigma2, 
    uni_types, types, pk, k_index, 
    gamma_1, gamma_1_rep_per_feature, gamma_2, 
    b, pi2, pi2_rep_per_feature, s2, Y, X, w_weight, v_weight,
    pi2_weighted_int, pi2_loglikeli, pi2_prop_n, MH_ind)

  ## Summary
  BETA <- rbind(rep(1, ncol(res$GAMMA_1)), res$GAMMA_1[k_index[-1],]) * 
    res$GAMMA_2 * res$B

  ## sum up duplicated features
  if(overlap_index){
    BETA_noInter <- BETA[-1, ]
    BETA_true <- sapply(1:nrow(U1), function(x) 
      apply(BETA_noInter[which(feature_dup_index==x), ,drop=FALSE], 2, sum))
    BETA <- rbind(BETA[1,], t(BETA_true))
  }
  random_sel <- sample(seq(1, nrow(BETA)), 5)
  beta_sub <- coda::as.mcmc(t(BETA[random_sel, ]))
  geweke_z <- coda::geweke.diag(beta_sub, frac1=0.1, frac2=0.5)$z
  convergence <- all(abs(geweke_z) < 1.96, na.rm=TRUE)

  totalIter <- keepIter
  # repeat if not convergence
  repeat{
    if(convergence | totalIter >= maxIter){
      break
    }
    cat("MCMC not converge. Running more ... \n")

    res <- MCMC_sb(seed, burnInIter=0, keepIter, 
      m1, P, printInt, N, 
      alpha1, alpha2, alpha_sigma, alpha_s, 
      beta1, beta2, beta_sigma, beta_s, pi1, sigma2, 
      uni_types, types, pk, k_index, 
      gamma_1, gamma_1_rep_per_feature, gamma_2, 
      b, pi2, pi2_rep_per_feature, s2, Y, X, w_weight, v_weight,
      pi2_weighted_int, pi2_loglikeli, pi2_prop_n, MH_ind)
    BETA <- rbind(rep(1, ncol(res$GAMMA_1)), res$GAMMA_1[k_index[-1],]) * 
      res$GAMMA_2 * res$B
    ## sum up duplicated features
    if(overlap_index){
      BETA_noInter <- BETA[-1, ]
      BETA_true <- sapply(1:nrow(U1), function(x) 
        apply(BETA_noInter[which(feature_dup_index==x), ,drop=FALSE], 2, sum))
      BETA <- rbind(BETA[1,], t(BETA_true))
    }
    random_sel <- sample(seq(1, nrow(BETA)), 5)
    beta_sub <- coda::as.mcmc(t(BETA[random_sel, ]))
    geweke_z <- coda::geweke.diag(beta_sub, frac1=0.1, frac2=0.5)$z
    convergence <- all(abs(geweke_z) < 1.96, na.rm=TRUE)
    totalIter <- totalIter + keepIter
  }
  res$totalIter <- totalIter

  ## Feature
  f_prob <- (apply(BETA!=0,1,mean))[-1] # exclude intercept
  f_null <- 1 - f_prob
  t <- seq(0,1,by=0.01)
  f_qvalue <- rep(NA,nrow(U1))
  for(j in 1:nrow(U1)){
    t_g <- t[t >= f_null[j]]
    f_qvalue[j] <- min(sapply(1:length(t_g),function(x)sum(f_null[f_null<=t_g[x]])/sum(f_null <= t_g[x])))
  }
  res$BETA <- BETA
  res$f_qvalue <- f_qvalue
  res$f_prob <- f_prob
  return(res)  
}
