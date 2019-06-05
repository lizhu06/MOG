#' MOG function for binary outcome.
#'
#' This function applies MOG to binary outcome
#'
#' @param X n by p feature matrix, where n is the number of samples, 
#' p is number of features
#' @param Y outcome vector of length n
#' @param U1 p by m1 matirx, denoting group membership of each feature
#' @param U2 m1 by m2 matirx, denoting group membership of each level-1 group
#' @param types numeric vector of length p, type index of each feature, allows different variance 
#' in priors of coefficient of different types 
#' @param lassoInit TRUE/FALSE, if TRUE, use lasso estimates as MCMC initials
#' @param weight_Tj TRUE/FALSE, if TRUE, the prior of gamma_jkl^0 will be weighted by T_j
#' @param weight_Dk TRUE/FALSE, if TRUE, the prior of gamma_kl^1 will be weighted by D_k
#' @param pi2_prop_n numeric number, specifying n_0 in the proposal distribution of pi_l^1
#' @param pi3_prop_n numeric number, specifying n_0 in the proposal distribution of pi_kl^0
#' @param MH_ind 1/0, if 1 MH is used to update pi_k^0; otherwise, Gibbs sampling is used
#' @param burnInIter number of iterations as burn-in period
#' @param keepIter number of iterations to be stored
#' @param maxIter maximum number of iterations
#' @param printInt print progress message every print_int iterations
#' @param seed integer, specifying the random seed
#' @param alpha1 numeric number, specifying alpha in the BETA distribution (prior of pi^2)
#' @param beta1 numeric number, specifying beta in the BETA distribution (prior of pi^2)
#' @param alpha2 numeric number, specifying alpha in the BETA distribution (prior of pi_l^1)
#' @param beta2 numeric number, specifying beta in the BETA distribution (prior of pi_l^1)
#' @param alpha3 numeric number, specifying alpha in the BETA distribution (prior of pi_kl^0)
#' @param beta3 numeric number, specifying beta in the BETA distribution (prior of pi_kl^0)
#' @return A list of MCMC output, together with feature selection probability 
#' and feature Bayesian qvalue
#' @export
MOG_binary <- function(X, Y, U1, U2, types=NULL, 
  lassoInit=FALSE,  weight_Tj=FALSE, weight_Dk=FALSE, 
  pi2_prop_n=10,  pi3_prop_n=10, MH_ind=1,
  burnInIter=1000, keepIter=3000, maxIter=10000, 
  printInt=100, seed=123, 
  alpha1=1, beta1=1, alpha2=1, beta2=1, alpha3=1, beta3=1){
  
  set.seed(seed)
  weight_Rj=FALSE
  weight_s2=FALSE
  # summarize the data
  N <- length(Y)
  if(is.null(types)){
    types <- rep(1, ncol(X))
  }
#  m2 <- ncol(U2)
#  m1 <- ncol(U1)
#  kl <- colSums(U2)
#  pk <- colSums(U1)

  # duplicate features if groups overlap
  U1_rowsum <- rowSums(U1)
  U2_rowsum <- rowSums(U2)
  if(sum(U1_rowsum==0)>1){
    stop("Feature must belongs to at least one group, make a singleton group for feature that does not belong sto any groups")
  }
  if(sum(U2_rowsum==0)>1){
    stop("Level-1 group must belongs to at least one level-2 group, make a singleton group for level-1 group that does not belong sto any level-2 groups")
  }

  overlap_index <- !(all(U1_rowsum==1) && all(U2_rowsum==1))
  if(overlap_index){
    cat("Data include overlapping groups. \n")
  }
  if(overlap_index){
    if(MH_ind==0){
      stop("When overlapping groups exist and Bernoulli needs to be corrected, MH_ind have to be set 1. \n")
    }
  }

  # create k_index, l_index, k_index_l
  feature_dup_index_U1 <- k_index_U1 <- NULL # corresponding to row index of U1 matrix
  for(j in 1:nrow(U1)){
    feature_dup_index_U1 <- c(feature_dup_index_U1, 
      rep(j, U1_rowsum[j]))
    k_index_U1 <- c(k_index_U1, which(U1[j,]==1))
  }

  # duplicate level-1 group (in the order of level-1 groups)
  feature_dup_index <- k_index <- l_index <- l_index_k <- NULL
  for(k in 1:nrow(U2)){
    k_counter <- length(unique(k_index)) 
    feature_dup_index <- c(feature_dup_index, 
      rep(feature_dup_index_U1[which(k_index_U1==k)], 
      times=U2_rowsum[k]))
    k_index <- c(k_index, rep(c((k_counter+1):(k_counter+U2_rowsum[k])),
      each=sum(k_index_U1==k)))
    l_index <- c(l_index, rep(which(U2[k, ]==1), each=sum(k_index_U1==k)))
    l_index_k <- c(l_index_k, which(U2[k, ]==1))
  }
  m1 <- length(unique(k_index))
  m2 <- length(unique(l_index))
  pk <- sapply(1:m1, function(k) sum(k_index==k))
  pl <- sapply(1:m2, function(l) sum(l_index==l))
  kl <- sapply(1:m2, function(l) sum(l_index_k==l))
  X_raw <- X
  X <- X_raw[, feature_dup_index]
  P <- ncol(X)
  types <- types[feature_dup_index]

  # add weight for duplicated features
  Rj <- sapply(1:length(feature_dup_index), function(x) 
    sum(feature_dup_index==feature_dup_index[x]))

  if(weight_s2){
    weight_for_s2 <- 1/Rj
  }else{
    weight_for_s2 <- rep(1, length(Rj))
  }

  if(weight_Rj){
    weight_for_Rj <- 1/Rj
  }else{
    weight_for_Rj <- rep(1, length(Rj))
  }

  if(weight_Dk){
    weight_for_Dk <- (1/U2_rowsum)
    weight_for_Dk <- rep(weight_for_Dk, times=U2_rowsum)
  }else{
    weight_for_Dk <- rep(1, length(U2_rowsum))
    weight_for_Dk <- rep(weight_for_Dk, times=U2_rowsum)
  }

  if(weight_Tj){
    weight_for_Tj <- (1/U1_rowsum)[feature_dup_index]
  }else{
    weight_for_Tj <- rep(1, length(Rj))
  }

  ##### hyper-parameters #####
  alpha_s<-0
  beta_s<-0
  alpha_sigma<-0
  beta_sigma<-0

  ## add intercept
  X <- cbind(rep(1, nrow(X)), X) # add intercept as the first column
  l_index <- c(m2+1, l_index) # add intercept
  k_index <- c(m1+1, k_index) # intercept have fake group index
  types <- c(max(types)+1, types) # intercept has fake type
  
  ##### set initial #####
  if(lassoInit == TRUE){
    ## Lasso 
    fitcv <- glmnet::cv.glmnet(X[,-1], Y, family="binomial", alpha=1, intercept=TRUE)
    coef <- coef(fitcv, s = "lambda.min")
    b <- coef[,1]
    gamma_1 <- sapply(1:m2, function(x) max((b!=0)[l_index==x]))
    gamma_2 <- sapply(1:m1, function(x) max((b!=0)[k_index==x]))
    gamma_3 <- b!=0
    pi1 <- max(sum(gamma_1)/length(gamma_1),0.1)
    pi2 <- sapply(1:m2, function(c) max(sum(gamma_2[l_index_k==c])/
      length(gamma_2[l_index_k==c]),0.1))
    pi3 <- sapply(1:m1, function(g) max(sum(gamma_3[k_index==g])/
      length(gamma_3[k_index==g]),0.1))
    pi2_rep_per_group <- pi2[l_index_k]
    pi3_rep_per_feature <- c(1, pi3[k_index[-1]])
    gamma_1_rep_per_feature <- c(1, gamma_1[l_index[-1]])
    gamma_2_rep_per_feature <- c(1, gamma_2[k_index[-1]])
    }else{
      b <- rnorm(P+1)
      pi1 <- 0.5
      pi2 <- rep(0.5, m2)
      pi2_rep_per_group <- pi2[l_index_k]
      pi3 <- rep(0.5,m1)
      pi3_rep_per_feature <- c(1, pi3[k_index[-1]])
      gamma_1 <- rbinom(m2,1,pi1)
      gamma_1_rep_per_feature <- c(1, gamma_1[l_index[-1]])
      gamma_2 <- rbinom(m1,1,pi2_rep_per_group)
      gamma_2_rep_per_feature <- c(1, gamma_2[k_index[-1]])
      gamma_3 <- rbinom(P+1,1,pi3_rep_per_feature)
    }
  sigma2 <- 1
  uni_types <- seq(1, max(types))  # including intercept
  s2 <- c(rep(1,length(uni_types)-1), 100)  # variance vector for coefficients of different types, last one is for intercept

  weight_s2_int <- weight_s2*1
  weight_Rj_int <- weight_Rj*1
  weight_Tj_int <- weight_Tj*1
  weight_Dk_int <- weight_Dk*1

  get_pi2_logl <- function(l, weight_Dk){
    sel_index <- which(l_index_k==l)
    if(weight_Dk){
      single_logl <- sum(gamma_2[sel_index]*log(weight_for_Dk[sel_index]*pi2[l])) + 
      sum((1-gamma_2[sel_index])*log(1-weight_for_Dk[sel_index]*pi2[l]))
    }else{
      single_logl <- sum(gamma_2[sel_index]*log(pi2[l])) + 
      sum((1-gamma_2[sel_index])*log(1-pi2[l]))
    }
    return(single_logl)
  }
  pi2_loglikeli <- sapply(1:m2, function(l) get_pi2_logl(l, weight_Dk))


  get_pi3_logl <- function(k, weight_Tj, weight_Rj){
    sel_index <- which(k_index==k)
    if(weight_Tj){
      single_logl <- sum(gamma_3[sel_index]*log(weight_for_Tj[sel_index]*pi3[k])) + 
      sum((1-gamma_3[sel_index])*log(1-weight_for_Tj[sel_index]*pi3[k]))
    }else if(weight_Rj){
      single_logl <- sum(gamma_3[sel_index]*log(weight_for_Rj[sel_index]*pi3[k])) + 
      sum((1-gamma_3[sel_index])*log(1-weight_for_Rj[sel_index]*pi3[k]))
    }else{
      single_logl <- sum(gamma_3[sel_index]*log(pi3[k])) + 
      sum((1-gamma_3[sel_index])*log(1-pi3[k]))
    }
    return(single_logl)
  }
  pi3_loglikeli <- sapply(1:m1, function(k) get_pi3_logl(k, weight_Tj, weight_Rj))


  #beta <- gamma_1_rep_per_feature*gamma_2_rep_per_feature*gamma_3*b

  res <- MCMC_mb(seed, burnInIter, keepIter, m2, m1, P, N, 
    alpha1, alpha2, alpha3, alpha_sigma, alpha_s, 
    beta1, beta2, beta3, beta_sigma, beta_s, pi1, sigma2, 
    uni_types, types, pl, pk, l_index, l_index_k, k_index, 
    gamma_1, gamma_1_rep_per_feature, gamma_2, 
    gamma_2_rep_per_feature, gamma_3, b, pi2, 
    pi2_rep_per_group, pi3, pi3_rep_per_feature, s2, Y, X, printInt, 
    weight_for_s2, weight_for_Rj, weight_for_Tj, weight_for_Dk,
    weight_s2_int, weight_Rj_int, weight_Tj_int, weight_Dk_int,
    pi2_loglikeli, pi3_loglikeli, pi2_prop_n, pi3_prop_n, MH_ind)

  ## Summary
  BETA_dup <- rbind(rep(1, ncol(res$GAMMA_1)),res$GAMMA_1[l_index[-1],])* 
   rbind(rep(1, ncol(res$GAMMA_1)),res$GAMMA_2[k_index[-1],]) * 
   res$GAMMA_3 * res$B

  ## sum up duplicated features
  BETA_noInter <- BETA_dup[-1, ]
  BETA_true <- sapply(1:nrow(U1), function(x) 
    apply(BETA_noInter[which(feature_dup_index==x), ,drop=FALSE], 2, sum))

  BETA <- rbind(BETA_dup[1,], t(BETA_true))
  rm(BETA_noInter)
  rm(BETA_true)

  ### check convergence
  random_sel <- sample(seq(1, nrow(BETA)), 5)
  beta_sub <- as.mcmc(t(BETA[random_sel, ]))
  geweke_z <- geweke.diag(beta_sub, frac1=0.1, frac2=0.5)$z
  convergence <- all(abs(geweke_z) < 1.96, na.rm=TRUE)

  totalIter <- keepIter

  ## repeat MCMC if it's not converged
  repeat{
    if(convergence | totalIter >= maxIter){
      break
    }
    cat("MCMC not converge. Running more ... \n")

    res <- MCMC_mb(seed, burnInIter=0, keepIter, m2, m1, P, N, 
      alpha1, alpha2, alpha3, alpha_sigma, alpha_s, 
      beta1, beta2, beta3, beta_sigma, beta_s, pi1, sigma2, 
      uni_types, types, pl, pk, l_index, l_index_k, k_index, 
      gamma_1, gamma_1_rep_per_feature, gamma_2, 
      gamma_2_rep_per_feature, gamma_3, b, pi2, 
      pi2_rep_per_group, pi3, pi3_rep_per_feature, s2, Y, X, printInt, 
      weight_for_s2, weight_for_Rj, weight_for_Tj, weight_for_Dk,
      weight_s2_int, weight_Rj_int, weight_Tj_int, weight_Dk_int,
      pi2_loglikeli, pi3_loglikeli, pi2_prop_n, pi3_prop_n, MH_ind)

    ## Summary
    BETA_dup <- rbind(rep(1, ncol(res$GAMMA_1)),res$GAMMA_1[l_index[-1],])* 
     rbind(rep(1, ncol(res$GAMMA_1)),res$GAMMA_2[k_index[-1],]) * 
     res$GAMMA_3 * res$B

    ## sum up duplicated features
    BETA_noInter <- BETA_dup[-1, ]
    BETA_true <- sapply(1:nrow(U1), function(x) 
      apply(BETA_noInter[which(feature_dup_index==x), ,drop=FALSE], 2, sum))

    BETA <- rbind(BETA_dup[1,], t(BETA_true))
    rm(BETA_noInter)
    rm(BETA_true)

    ### check convergence
    random_sel <- sample(seq(1, nrow(BETA)), 5)
    beta_sub <- as.mcmc(t(BETA[random_sel, ]))
    geweke_z <- geweke.diag(beta_sub, frac1=0.1, frac2=0.5)$z
    convergence <- all(abs(geweke_z) < 1.96, na.rm=TRUE)

    totalIter <- keepIter + totalIter

  }

  ## Feature
  f_ind <- (apply(BETA!=0,1,mean))[-1] # exclude intercept
  f_null <- 1-f_ind
  t <- seq(0,1,by=0.01)
  f_qvalue <- rep(NA,nrow(U1))
  for(j in 1:nrow(U1)){
    t_g <- t[t >= f_null[j]]
    f_qvalue[j] <- min(sapply(1:length(t_g),function(x)sum(f_null[f_null<=t_g[x]])/sum(f_null<=t_g[x])))
  }
  res$feature_dup_index <- feature_dup_index # index in colnames(X) and rownames(U1)
  res$k_index <- k_index
  res$l_index <- l_index

  res$BETA <- BETA
  res$BETA_dup <- BETA_dup
  res$f_qvalue <- f_qvalue
  res$f_prob <- f_ind  

  res$pi2_loglikeli <- pi2_loglikeli
  res$pi3_loglikeli <- pi3_loglikeli
  res$weight_for_Rj <- weight_for_Rj
  res$weight_for_Dk <- weight_for_Dk
  res$weight_for_Tj <- weight_for_Tj
  return(res)  
}
