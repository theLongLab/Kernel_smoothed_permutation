library(utilities)
library(MASS)
library(moments)

chi_square_significant_genotype <- read.csv("chi_square_significant_genotype.csv", header = F)
## determine the SNP you want to do permutaion
SNP_index <- 1 # random, can choose your own
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

## permutation, start with data_standard (removed records with genotype=9)
standard <- floor(log10(1/pvalue_standard))+1
sim_time <- 1.1*10^standard

chi.square.perm <- function(data_standard) {
  ## permutate phenotype: not include missing genotype records
  y_perm <- sample(data_standard$pheno, dim(data_standard)[1], replace = FALSE)
  X_tmp <- chisq.test(data_standard$geno, y_perm)$statistic
  names(X_tmp) <- NULL
  return(X_tmp)
}

X_perm <- replicate(sim_time, chi.square.perm(data_standard))
write.table(X_perm,paste("X_perm_", chr_target, "_", loc_target, ".txt", sep = ""), sep="\t", row.names=F, quote=F, col.names=F)




