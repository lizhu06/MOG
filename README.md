# MOG

R package of Bayesian indicator variable selection incorporating multi-layer overlapping groups

## Required Package

glmnet, coda, Rcpp, RcppEigen

### Installing

In R console

```
library(devtools)
install_github("lizhu06/MOG")
```

## Running a simulation
```
library(MOG)

data <- geneSimu1(seed=1)
res <- SOG_continuous(data$X, data$Y, data$U1, 
                      lassoInit=FALSE, BernoulliWeighted=FALSE, 
                      pi2_prop_n=10, MH_ind=1,
                      burnInIter=1000, keepIter=3000, maxIter=5000, 
                      printInt=1000, seed=123, 
                      PreEstBeta=FALSE, alpha2=1, beta2=1)
feature_FDR <-sum(data$true_beta==0 & res$f_qvalue<=0.1)/sum(res$f_qvalue<=0.1)
feature_FOR <-sum(data$true_beta!=0 & res$f_qvalue>0.1)/sum(res$f_qvalue>0.1)
AUC <-pROC::auc(pROC::roc(response=factor(data$true_beta!=0), 
                          predictor=res$f_prob))
beta_median <- apply(res$BETA[,-seq(1,1000)],1,median)
MSE <- mean((data$X %*% beta_median - data$Y)^2)

```

## Authors

* **Li Zhu** - liz86@pitt.edu


