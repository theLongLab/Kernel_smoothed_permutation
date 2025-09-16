

############################ GWAS KSP ########################################
setwd("your_folder_path")
source("kernelSmoothedPerm.R")

X_perm <- read.table("snp1_gwas_perm.txt", sep = "\t")
X_perm <- X_perm$V1
print(summary(X_perm))
print(length(X_perm))

X_standard <- 5.26205200943128
pvalue_threshold=1e-07

results <- kernelSmoothedPerm(X_standard, X_perm, pvalue_threshold, sim_time = 10)
print(results)

# $optimal_lambda
# [1] "none"
# 
# $pvalue_naive
# [1] 8.810573e-08
# 
# $pvalue_kernel
# [1] 1.525465e-07

write.csv(results,file = "snp1_KSP_results.csv")



