library(utilities)
library(MASS)
library(moments)

chi_square_significant_genotype <- read.csv("chi_square_significant_genotype.csv", header = F)
## determine the SNP you want to do permutaion
SNP_index <- 1
number_of_indiv <- 4798 # 1860 case + 2938 control

SNP_target <- chi_square_significant_genotype[SNP_index, ]
chr_target <- unlist(SNP_target, use.names = FALSE)[1]
loc_target <- unlist(SNP_target, use.names = FALSE)[2]
genotype_target <- unlist(SNP_target, use.names = FALSE)[3:(number_of_indiv+2)]

pheno <- read.table(file='RA.tsv')
pheno$V3 <- ifelse(pheno$V2==1,"control","case")
y.b <- pheno$V3

## standard chi-squared test
data_standard <- data.frame(genotype_target, y.b)
colnames(data_standard) <- c("geno","pheno")
data_standard <- data_standard[which(data_standard$geno != 9),]

pvalue_standard <- chisq.test(data_standard$geno,data_standard$pheno)$p.value
X_standard <- chisq.test(data_standard$geno,data_standard$pheno)$statistic
names(X_standard) <- NULL

print(pvalue_standard)
print(X_standard)

## sim_time
standard <- floor(log10(1/pvalue_standard))+1
sim_time <- 10^standard
## repeat_time_in: for each setting, we repeat for 300 times and take the median
repeat_time_in <- 300
## band_coef: bandwith coefficients we want to test
band_coef <- seq(0.5,15,0.5)
#length(band_coef) # 30

## read permutated test statistics
X_perm <- read.table(paste("X_perm_", chr_target, "_", loc_target, ".txt", sep = ""),sep = "\t")
X_perm <- X_perm$V1

## KDE
## repeat for 300 times and take median
for (i in 1:repeat_time_in){
  kernel_result <- data.frame(matrix(0, repeat_time_in, 5+length(band_coef)*2+4))
  colnames(kernel_result) <- c("chr",
                               "loc",
                               "pvalue_standard", 
                               "X_standard",
                               "pvalue_naive", 
                               paste("pvalue_KDE_full_",band_coef,sep=""),
                               paste("pvalue_KDE_sub_",band_coef,sep=""),
                               "skewness_full",
                               "kurtosis_full",
                               "skewness_sub",
                               "kurtosis_sub")
  
  kernel_result[i,1] <- chr_target
  kernel_result[i,2] <- loc_target
  kernel_result[i,3] <- pvalue_standard
  kernel_result[i,4] <- X_standard
  
  ## prepare data
  X_perm_each <- sample(X_perm, sim_time, replace = F)
  X_perm_each_sub <- sample(X_perm_each,sim_time/10,replace = F) # check 10% sub-samples
  
  ##### step 2: naive permutation
  pvalue_naive <- sum(X_perm_each >= X_standard)/sim_time
  kernel_result[i,5] <- pvalue_naive
  
  ##### step 3: perform KDE
  ## 1: Original transformation
  X_perm_each_transform <- X_perm_each
  X_standard_transform <- X_standard
  X_perm_each_sub_transform <- X_perm_each_sub
  
  for(j in 1:length(band_coef)){
    MY_KDE1 <- KDE(X_perm_each_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_transform), to.environment = T) 
    kernel_result[i,5+j] <- pkde(X_standard_transform,lower.tail = F) 
    MY_KDE2 <- KDE(X_perm_each_sub_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_sub_transform), to.environment = T) 
    kernel_result[i,5+length(band_coef)+j] <- pkde(X_standard_transform,lower.tail = F)
  }
  
  kernel_result[i,5+length(band_coef)*2+1] <- skewness(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+2] <- kurtosis(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+3] <- skewness(X_perm_each_sub_transform)
  kernel_result[i,5+length(band_coef)*2+4] <- kurtosis(X_perm_each_sub_transform)
  
  write.csv(kernel_result, "chi_square_KDE_original.csv", row.names=FALSE)
  
  ## 2: Squared root transformation
  X_perm_each_transform <- sqrt(X_perm_each)
  X_standard_transform <- sqrt(X_standard)
  X_perm_each_sub_transform <- sqrt(X_perm_each_sub)
  
  for(j in 1:length(band_coef)){
    MY_KDE1 <- KDE(X_perm_each_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_transform), to.environment = T) 
    kernel_result[i,5+j] <- pkde(X_standard_transform,lower.tail = F) 
    MY_KDE2 <- KDE(X_perm_each_sub_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_sub_transform), to.environment = T) 
    kernel_result[i,5+length(band_coef)+j] <- pkde(X_standard_transform,lower.tail = F)
  }
  
  kernel_result[i,5+length(band_coef)*2+1] <- skewness(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+2] <- kurtosis(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+3] <- skewness(X_perm_each_sub_transform)
  kernel_result[i,5+length(band_coef)*2+4] <- kurtosis(X_perm_each_sub_transform)
  
  write.csv(kernel_result, "chi_square_KDE_square_root.csv", row.names=FALSE)
  
  ## 3: Cubic root transformation
  X_perm_each_transform <- X_perm_each^(1/3)
  X_standard_transform <- X_standard^(1/3)
  X_perm_each_sub_transform <- X_perm_each_sub^(1/3)
  
  for(j in 1:length(band_coef)){
    MY_KDE1 <- KDE(X_perm_each_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_transform), to.environment = T) 
    kernel_result[i,5+j] <- pkde(X_standard_transform,lower.tail = F) 
    MY_KDE2 <- KDE(X_perm_each_sub_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_sub_transform), to.environment = T) 
    kernel_result[i,5+length(band_coef)+j] <- pkde(X_standard_transform,lower.tail = F)
  }
  
  kernel_result[i,5+length(band_coef)*2+1] <- skewness(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+2] <- kurtosis(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+3] <- skewness(X_perm_each_sub_transform)
  kernel_result[i,5+length(band_coef)*2+4] <- kurtosis(X_perm_each_sub_transform)
  
  write.csv(kernel_result, "chi_square_KDE_cubic_root.csv", row.names=FALSE)
  
  ## 4: Log transformation
  X_perm_each_transform <- log(X_perm_each,base = exp(1)) 
  X_standard_transform <- log(X_standard,base = exp(1)) 
  X_perm_each_sub_transform <- log(X_perm_each_sub,base = exp(1))
  
  for(j in 1:length(band_coef)){
    MY_KDE1 <- KDE(X_perm_each_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_transform), to.environment = T) 
    kernel_result[i,5+j] <- pkde(X_standard_transform,lower.tail = F) 
    MY_KDE2 <- KDE(X_perm_each_sub_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_sub_transform), to.environment = T) 
    kernel_result[i,5+length(band_coef)+j] <- pkde(X_standard_transform,lower.tail = F)
  }
  
  kernel_result[i,5+length(band_coef)*2+1] <- skewness(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+2] <- kurtosis(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+3] <- skewness(X_perm_each_sub_transform)
  kernel_result[i,5+length(band_coef)*2+4] <- kurtosis(X_perm_each_sub_transform)
  
  write.csv(kernel_result, "chi_square_KDE_log.csv", row.names=FALSE)
  
  ## 5: Box-Cox transformation
  ## find optimal power
  cox_X <- boxcox(lm(X_perm_each ~ 1),plotit = F)
  lambda_X <- cox_X$x[which.max(cox_X$y)]
  
  ## from here, we add 1 more column on kernel_result to represent optimal lambda
  if(lambda_X == 0){
    lambda_X_optinal <- lambda_X
    ## Box-Cox transformation = Log transformation
    kernel_result$lambda_optinal <- lambda_X_optinal
    write.csv(kernel_result, "chi_square_KDE_Box.csv", row.names=FALSE)
  } else {
    lambda_X_optinal <- lambda_X
    ## Box-Cox transformation
    X_perm_each_transform <- (X_perm_each^lambda_X_optinal-1)/lambda_X_optinal
    X_standard_transform <- (X_standard^lambda_X_optinal-1)/lambda_X_optinal
    X_perm_each_sub_transform <- (X_perm_each_sub^lambda_X_optinal-1)/lambda_X_optinal
    
    for(j in 1:length(band_coef)){
      MY_KDE1 <- KDE(X_perm_each_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_transform), to.environment = T) 
      kernel_result[i,5+j] <- pkde(X_standard_transform,lower.tail = F) 
      MY_KDE2 <- KDE(X_perm_each_sub_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_sub_transform), to.environment = T) 
      kernel_result[i,5+length(band_coef)+j] <- pkde(X_standard_transform,lower.tail = F)
    }
    
    kernel_result[i,5+length(band_coef)*2+1] <- skewness(X_perm_each_transform)
    kernel_result[i,5+length(band_coef)*2+2] <- kurtosis(X_perm_each_transform)
    kernel_result[i,5+length(band_coef)*2+3] <- skewness(X_perm_each_sub_transform)
    kernel_result[i,5+length(band_coef)*2+4] <- kurtosis(X_perm_each_sub_transform)
    
    kernel_result$lambda_optinal <- lambda_X_optinal
    write.csv(kernel_result, "chi_square_KDE_Box.csv", row.names=FALSE)
  }
  
  ## 6: Skewness-driven Box-Cox transformation
  ## find optimal Skewness-driven power: as close as 0
  lambda_X_test <- lambda_X-seq(-0.2,0.2,0.01)
  lambda_X_test <- lambda_X_test[lambda_X_test != 0]
  
  Skewness_X_test <- numeric(0)
  for (k in 1:length(lambda_X_test)) {
    X_test <- (X_perm_each^lambda_X_test[k]-1)/lambda_X_test[k]
    Skewness_X_test[k] <- skewness(X_test)
  }
  lambda_X_optinal <- lambda_X_test[which.min(abs(Skewness_X_test))]
  
  ## Box-Cox transformation
  X_perm_each_transform <- (X_perm_each^lambda_X_optinal-1)/lambda_X_optinal
  X_standard_transform <- (X_standard^lambda_X_optinal-1)/lambda_X_optinal
  X_perm_each_sub_transform <- (X_perm_each_sub^lambda_X_optinal-1)/lambda_X_optinal
  
  for(j in 1:length(band_coef)){
    MY_KDE1 <- KDE(X_perm_each_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_transform), to.environment = T) 
    kernel_result[i,5+j] <- pkde(X_standard_transform,lower.tail = F) 
    MY_KDE2 <- KDE(X_perm_each_sub_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_sub_transform), to.environment = T) 
    kernel_result[i,5+length(band_coef)+j] <- pkde(X_standard_transform,lower.tail = F)
  }
  
  kernel_result[i,5+length(band_coef)*2+1] <- skewness(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+2] <- kurtosis(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+3] <- skewness(X_perm_each_sub_transform)
  kernel_result[i,5+length(band_coef)*2+4] <- kurtosis(X_perm_each_sub_transform)
  
  kernel_result$lambda_optinal <- lambda_X_optinal
  write.csv(kernel_result, "chi_square_KDE_Box_Skew.csv", row.names=FALSE)
  
  ## 7: Kurtosis-driven Box-Cox transformation
  ## find optimal kurtosis-driven power: as close as 3
  lambda_X_test <- lambda_X-seq(-0.2,0.2,0.01)
  lambda_X_test <- lambda_X_test[lambda_X_test != 0]
  
  kurtosis_X_test <- numeric(0)
  for (k in 1:length(lambda_X_test)) {
    X_test <- (X_perm_each^lambda_X_test[k]-1)/lambda_X_test[k]
    kurtosis_X_test[k] <- kurtosis(X_test)
  }
  lambda_X_optinal <- lambda_X_test[which.min(abs(kurtosis_X_test-3))]
  
  ## Box-Cox transformation
  X_perm_each_transform <- (X_perm_each^lambda_X_optinal-1)/lambda_X_optinal
  X_standard_transform <- (X_standard^lambda_X_optinal-1)/lambda_X_optinal
  X_perm_each_sub_transform <- (X_perm_each_sub^lambda_X_optinal-1)/lambda_X_optinal
  
  for(j in 1:length(band_coef)){
    MY_KDE1 <- KDE(X_perm_each_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_transform), to.environment = T) 
    kernel_result[i,5+j] <- pkde(X_standard_transform,lower.tail = F) 
    MY_KDE2 <- KDE(X_perm_each_sub_transform, bandwidth = band_coef[j]*bw.nrd0(X_perm_each_sub_transform), to.environment = T) 
    kernel_result[i,5+length(band_coef)+j] <- pkde(X_standard_transform,lower.tail = F)
  }
  
  kernel_result[i,5+length(band_coef)*2+1] <- skewness(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+2] <- kurtosis(X_perm_each_transform)
  kernel_result[i,5+length(band_coef)*2+3] <- skewness(X_perm_each_sub_transform)
  kernel_result[i,5+length(band_coef)*2+4] <- kurtosis(X_perm_each_sub_transform)
  
  kernel_result$lambda_optinal <- lambda_X_optinal
  write.csv(kernel_result, "chi_square_KDE_Box_Kurt.csv", row.names=FALSE)
  
}




