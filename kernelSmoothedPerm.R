####################################################################
library(utilities)
library(MASS)
library(moments)

kernelSmoothedPerm <- function(X_standard, X, pvalue_threshold){
  ### 0. P-value estimation for naive permutation
  pvalue_naive = sum(X >= X_standard)/length(X)
  
  ### 1. check kurtosis of permuted test statistics
  if(abs(kurtosis(X)-3) <= 0.01){
    lambda_X_optinal = "none"
    X_trans = X
    X_standard_trans = X_standard  ## go to 4
  } else {
    ### 2. Box-Cox transformation: data must be positive (no negative values) and continuous
    ###    If there is negative value, then data should be normalized.
    if(sum(X <= 0) > 0){
      X = (X-min(X))/(max(X)-min(X))
      X[which(X == 0)] = 1e-05
    }
    cox_X = boxcox(lm(X ~ 1),plotit = F)
    lambda_X = cox_X$x[which.max(cox_X$y)]
    
    ### 3. find the optimal power and transform the data
    ## 3.1 if optimal power = 0, then perform log transformation
    if(lambda_X == 0){
      lambda_X_optinal = 0
      X_trans = log(X,base = exp(1))
      X_standard_trans = log(X_standard,base = exp(1))  ## go to 4
    } else {
      ## 3.2 find the optimal kurtosis-driven power and transform the data
      lambda_X_test = lambda_X-seq(-0.2,0.2,0.01)
      lambda_X_test = lambda_X_test[lambda_X_test != 0]
      kurtosis_X_test = numeric(0)
      for (i in 1:length(lambda_X_test)) {
        X_trans_test = (X^lambda_X_test[i]-1)/lambda_X_test[i]
        kurtosis_X_test[i] = kurtosis(X_trans_test)
      }
      lambda_X_optinal = lambda_X_test[which.min(abs(kurtosis_X_test-3))]
      X_trans = (X^lambda_X_optinal-1)/lambda_X_optinal
      X_standard_trans = (X_standard^lambda_X_optinal-1)/lambda_X_optinal ## go to 4
    }
  }
  
  ### 4. apply KDE based on P-value accuracy threshold
  if(pvalue_threshold == 1e-07){
    pvalue_test = numeric(0)
    for (j in 1:50){
      X_trans_test = sample(X_trans,1e+07,replace = F)
      MY_KDE = KDE(X_trans_test,bandwidth = 5*bw.nrd0(X_trans_test),to.environment = T)
      pvalue_test[j] = pkde(X_standard_trans,lower.tail = F)
    }
    pvalue_kernel = median(pvalue_test)
  }  
  
  if(pvalue_threshold == 1e-08){
    pvalue_test = numeric(0)
    for (j in 1:50){
      X_trans_test = sample(X_trans,1e+08,replace = F)
      MY_KDE = KDE(X_trans_test,bandwidth = 7*bw.nrd0(X_trans_test),to.environment = T)
      pvalue_test[j] = pkde(X_standard_trans,lower.tail = F)
    }
    pvalue_kernel = median(pvalue_test)
  }
  
  newlist = list(optimal_lambda = lambda_X_optinal, pvalue_naive = pvalue_naive, pvalue_kernel = pvalue_kernel)
  return(newlist)
}
