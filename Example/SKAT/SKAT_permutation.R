library(SKAT)
library(utilities)
## determine the gene you want to do permutaion
#gene_name = "gene6281"
#setwd(paste("/work/long_lab/jbian/Kernel/WG_SKAT/WG_SKAT_RA/",gene_name,"/",sep=""))

## X: covariate
number_of_indiv = 4798 # 1860 case + 2938 control
X = rep(1,number_of_indiv)
X = as.matrix(X)
## y.b: phenotype
pheno = read.table(file="RA.tsv")
pheno$V2 = ifelse(pheno$V2==1,0,1)
table(pheno$V2)
# 0    1 
# 2938 1860 
y.b = pheno$V2
## Z: geno matrix 
data_input = read.csv("SKAT_Z.csv", header = F)
data_input = data_input[,3:(number_of_indiv+2)]
data_input = as.matrix(data_input)
data_input = t(data_input)
Z = data_input

## SKAT model
obj <- SKAT_Null_Model(y.b ~ X, out_type="D")
X_standard <- SKAT(Z, obj, kernel = "linear.weighted")$Q 
pvalue_standard = SKAT(Z, obj, kernel = "linear.weighted")$p.value 
print(X_standard)
print(pvalue_standard)

## permutation
standard = floor(log10(1/pvalue_standard))+1
sim_time = 1.1*10^standard

SKAT.perm <- function(y.b, X, Z) {
  ## permutate phenotype
  y_perm = sample(y.b,number_of_indiv,replace = FALSE)
  obj <- SKAT_Null_Model(y_perm ~ X, out_type="D") 
  X_tmp = SKAT(Z, obj, kernel = "linear.weighted")$Q
  return(X_tmp)
}

X_perm = replicate(sim_time, SKAT.perm(y.b, X, Z))
write.table(X_perm, "X_perm.txt", sep = "\t", row.names = F, col.names = F)
## Note!!!: in reality generate all permutated test-statistics takes a lot of time, 
## especially for pvalue level at e-08, so we need to generate small samples each time 
## and combine the samples as final samples




