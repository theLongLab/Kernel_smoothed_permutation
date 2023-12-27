####################################################################
library(utilities)
library(MASS)
library(moments)

kernelSmoothedPerm <- function(a_standard,a,pvalue_threshold){
  ### 0. P-value estimation for naive permutation
  pvalue_naive = sum(a >= a_standard)/length(a)
  
  ### 1. check kurtosis of permuted test statistics
  if(abs(kurtosis(a)-3) <= 0.01){
    lambda_a_optinal = "none"
    a_trans = a
    a_standard_trans = a_standard  ## go to 4
  } else {
    ### 2. Box-Cox transformation: data must be positive (no negative values) and continuous
    ###    If there is negative value, then data should be normalized.
    if(sum(a <= 0) > 0){
      a = (a-min(a))/(max(a)-min(a))
      a[which(a == 0)] = 1e-05
    }
    cox_a = boxcox(lm(a ~ 1),plotit = F)
    lambda_a = cox_a$x[which.max(cox_a$y)]
    
    if(lambda_a == 0){
      lambda_a_optinal = 0
      a_trans = log(a,base = exp(1))
      a_standard_trans = log(a_standard,base = exp(1))  ## go to 4
    } else {
      ### 3. find the optimal kurtosis-driven power and transform the data
      lambda_a_test = lambda_a-seq(-0.2,0.2,0.01)
      lambda_a_test = lambda_a_test[lambda_a_test != 0]
      kurtosis_a_test = numeric(0)
      for (i in 1:length(lambda_a_test)) {
        a_trans_test = (a^lambda_a_test[i]-1)/lambda_a_test[i]
        kurtosis_a_test[i] = kurtosis(a_trans_test)
      }
      lambda_a_optinal = lambda_a_test[which.min(abs(kurtosis_a_test-3))]
      a_trans = (a^lambda_a_optinal-1)/lambda_a_optinal
      a_standard_trans = (a_standard^lambda_a_optinal-1)/lambda_a_optinal ## go to 4
    }
  }
  
  ### 4. apply KDE based on P-value accuracy threshold
  if(pvalue_threshold == 1e-07){
    pvalue_test = numeric(0)
    for (j in 1:50){
      a_trans_test = sample(a_trans,1e+07,replace = F)
      MY_KDE = KDE(a_trans_test,bandwidth = 5*bw.nrd0(a_trans_test),to.environment = T)
      pvalue_test[j] = pkde(a_standard_trans,lower.tail = F)
    }
    pvalue_kernel = median(pvalue_test)
  }  
  
  if(pvalue_threshold == 1e-08){
    pvalue_test = numeric(0)
    for (j in 1:50){
      a_trans_test = sample(a_trans,1e+08,replace = F)
      MY_KDE = KDE(a_trans_test,bandwidth = 7*bw.nrd0(a_trans_test),to.environment = T)
      pvalue_test[j] = pkde(a_standard_trans,lower.tail = F)
    }
    pvalue_kernel = median(pvalue_test)
  }
  
  newlist = list(optimal_lambda = lambda_a_optinal, pvalue_naive = pvalue_naive, pvalue_kernel = pvalue_kernel)
  return(newlist)
}







