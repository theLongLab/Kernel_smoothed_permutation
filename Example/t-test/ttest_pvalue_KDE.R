# install.packages("remotes")   # lightweight helper
# remotes::install_github("ben-oneill/utilities")
library(utilities)

##### problem setting: two-sample t-test 
## pvalue_level, sim_time: 
## if target p-value level is e-07, then set pvalue_level as 7
## if target p-value level is e-08, then set pvalue_level as 8
pvalue_level = 7
sim_time = 10^pvalue_level
## repeat_time_in: for each setting, we repeat for 50 times and take the median
## repeat_time_out: we run analysis for 500 settings
repeat_time_in = 50 
repeat_time_out = 300
## band_coef: bandwith coefficients we want to test
band_coef = seq(0.5,15,0.5)
#length(band_coef) # 30
## kernel_result
kernel_result = data.frame(matrix(0, repeat_time_out, 8+length(band_coef)*2))
colnames(kernel_result) = c("sample size", 
                            "mean a", 
                            "mean b", 
                            "sd a", 
                            "sd b",
                            "pvalue_standard", 
                            "X_standard",
                            "pvalue_naive", 
                            paste("pvalue_KDE_full_",band_coef,sep=""),
                            paste("pvalue_KDE_sub_",band_coef,sep=""))

##### t.test.perm: function of doing permutation for t-test
## x : combined two random samples
t.test.perm <- function(x, na, nb) {
  tmp_sample = sample(x, (na+nb),replace = FALSE)
  tmp_a = tmp_sample[1:na]
  tmp_b = tmp_sample[(1+na):(na+nb)]
  sp_2 = ((na-1)*var(tmp_a)+(nb-1)*var(tmp_b))/(na+nb-2)
  X_tmp = (mean(tmp_a)-mean(tmp_b))/(sqrt(sp_2)*sqrt(1/na+1/nb))
  return(X_tmp)
}

##### 
for (i in 1:repeat_time_out){
  ##### step 1: standard test
  na = sample(c(100,150,200,250,500),1,replace = F)
  nb = na
  mean_a = runif(1,0,1)
  mean_b = 0
  sd_a = sample(c(1,2,3,4,5),1,replace = F)
  sd_b = sd_a
  
  ## pvalue_standard
  repeat{
    random_a = rnorm(na,mean_a,sd_a)
    random_b = rnorm(nb,mean_b,sd_b) 
    standard_test = t.test(random_a, random_b, alternative = "two.sided",
                           mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
    pvalue_standard = standard_test$p.value 
    standard = floor(log10(1/pvalue_standard))+1
    if(standard == pvalue_level){
      break
    }
  }
  
  kernel_result[i,1] = na
  kernel_result[i,2] = mean_a
  kernel_result[i,3] = mean_b
  kernel_result[i,4] = sd_a
  kernel_result[i,5] = sd_b
  kernel_result[i,6] = pvalue_standard
  
  ## X_standard
  sp_2 = ((na-1)*var(random_a)+(nb-1)*var(random_b))/(na+nb-2)
  X_standard = (mean(random_a)-mean(random_b))/(sqrt(sp_2)*sqrt(1/na+1/nb))
  kernel_result[i,7] = X_standard
  
  ##### step 2: generate permutated test statistics
  X_perm = replicate(1.1*sim_time, t.test.perm(x = c(random_a,random_b), na = na, nb = nb))
  write.table(X_perm, paste("X_perm", i, ".txt", sep=""), sep = "\t", row.names = F, col.names = F)
  
  ##### repeat step 3 and 4 for 50 times and take median
  ## we also keep track of original results for each setting tested
  kernel_result_tmp = data.frame(matrix(0, repeat_time, 8+length(band_coef)*2))
  colnames(kernel_result_tmp) = colnames(kernel_result)
  kernel_result_tmp[,1] = na
  kernel_result_tmp[,2] = mean_a
  kernel_result_tmp[,3] = mean_b
  kernel_result_tmp[,4] = sd_a
  kernel_result_tmp[,5] = sd_b
  kernel_result_tmp[,6] = pvalue_standard
  kernel_result_tmp[,7] = X_standard
  
  for (j in 1:repeat_time_in){
    ## prepare data
    X_perm_each = sample(X_perm, sim_time, replace = F)
    X_perm_each_sub = sample(X_perm_each,sim_time/10,replace = F) # check 10% sub-samples
    ##### step 3: naive permutation
    X_all = c(X_standard, X_perm_each)
    X_all = sort(X_all)
    if(mean_a > mean_b){
      pvalue_naive = 1-which(X_all==X_standard)/(sim_time+1)
    }else{
      pvalue_naive = which(X_all==X_standard)/(sim_time+1)
    }
    kernel_result_tmp[j,8] = pvalue_naive
    
    ##### step 4: KDE
    for(k in 1:length(band_coef)){
      MY_KDE1 = KDE(X_perm_each, bandwidth = band_coef[k]*bw.nrd0(X_perm_each), to.environment = T) 
      kernel_result_tmp[j,8+k] = pkde(X_standard,lower.tail = F) 
      MY_KDE2 = KDE(X_perm_each_sub, bandwidth = band_coef[k]*bw.nrd0(X_perm_each_sub), to.environment = T) 
      kernel_result_tmp[j,8+length(band_coef)+k] = pkde(X_standard,lower.tail = F)
    }
  }
  write.csv(kernel_result_tmp, paste("KDE_raw", i, ".csv", sep=""), row.names=FALSE)
  
  ## update kernel_result using median
  kernel_result[i, 9:dim(kernel_result_tmp)[2]] = as.vector(sapply(kernel_result_tmp, median))[9:dim(kernel_result_tmp)[2]]
}

write.table(kernel_result,"KDE_ttest_e07.csv",sep=",",row.names=FALSE)

