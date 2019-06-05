# MOG

R package of Bayesian indicator variable selection incorporating multi-layer overlapping groups

## Required Package

glmnet, Rcpp, RcppEigen

### Installing

In R console

```
library(devtools)
install_github("lizhu06/MOG")
```


## Running the simulations

### Simulation 1

```
library(MOG)

data <- geneSimu1(seed=1)
res <- SOG_continuous(data$X, data$Y, data$U1, totalIter=2000,
                      burnInIter=1000, seed=123, printInt=100)
feature_FDR <-sum(data$true_beta==0 & res$f_qvalue<=0.1)/sum(res$f_qvalue<=0.1)
feature_FOR <-sum(data$true_beta!=0 & res$f_qvalue>0.1)/sum(res$f_qvalue>0.1)
AUC <-pROC::auc(pROC::roc(response=factor(data$true_beta!=0), 
                          predictor=res$f_prob))
beta_median <- apply(res$BETA[,-seq(1,1000)],1,median)
MSE <- mean((data$X %*% beta_median - data$Y)^2)

```

### Simulation 2

```
data <- geneSimu2(seed=1)
res <- SOG_continuous(data$X, data$Y, data$U1, totalIter=2000,
                      burnInIter=1000, seed=123, printInt=100)
feature_FDR <-sum(data$true_beta==0 & res$f_qvalue<=0.1)/sum(res$f_qvalue<=0.1)
feature_FOR <-sum(data$true_beta!=0 & res$f_qvalue>0.1)/sum(res$f_qvalue>0.1)
AUC <-pROC::auc(pROC::roc(response=factor(data$true_beta!=0), 
                          predictor=res$f_prob))
beta_median <- apply(res$BETA[,-seq(1,1000)],1,median)
MSE <- mean((data$X %*% beta_median - data$Y)^2)

```

### Simulation 3

```
data <- geneSimu3(U=0.5, seed=1)
res <- MOG_continuous(data$X, data$Y, data$U1, data$U2, totalIter=2000,
                      burnInIter=1000, seed=123, printInt=100)
feature_FDR <-sum(data$true_beta==0 & res$f_qvalue<=0.1)/sum(res$f_qvalue<=0.1)
feature_FOR <-sum(data$true_beta!=0 & res$f_qvalue>0.1)/sum(res$f_qvalue>0.1)
AUC <-pROC::auc(pROC::roc(response=factor(data$true_beta!=0), 
                          predictor=res$f_prob))
beta_median <- apply(res$BETA[,-seq(1,1000)],1,median)
MSE <- mean((data$X %*% beta_median - data$Y)^2)

```
## Authors

* **Li Zhu** - liz86@pitt.edu


