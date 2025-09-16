
## ===== Read: One SNP genotype + binary phenotype =====
## Genotype file: geno.csv (column containing only one SNP)
## Phenotype file: pheno.txt (two columns: FID, PHENO (0/1))

# --- read genotype (keep a single SNP column) ---
G=read.csv("geno.csv")
rownames(G)=G[,1]
G=G[,-c(1:5)]
G <- t(as.matrix(G))

g <- as.numeric(G[, 1])            # <== Modify to your single SNP column index or column name
snp_name <- colnames(G)[1]         # <== Same as above

# --- read phenotype ---
ph <- read.table("pheno.txt", header = FALSE, stringsAsFactors = FALSE)
rownames(ph) <- ph[,1]
y <- as.vector(ph[,2])
stopifnot(length(g) == length(y))
stopifnot(all(y %in% c(0,1)))

n <- length(y)

## ========= Parameter settings =========
B_total      <- 1.1*1e7            # Number of replacements
batch_size   <- 10000        # Batch size for each disk write (per core)
out_dir      <- "perm_out"     # Output directory (one file per core)
use_gzip     <- FALSE           # Whether to write gzip compressed csv
ncores       <- min(60L, parallel::detectCores() - 1L)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========= Single fitting function: returns c(Z, OR, P_add) =========
fit_once <- function(yv, gv) {
  # Quality threshold: at least two genotypes and sufficient samples
  if (length(unique(gv)) < 2 || length(yv) < 10) return(c(NA_real_, NA_real_, NA_real_))
  # Fitting, capturing separation/non-convergence issues
  res <- tryCatch({
    fit <- suppressWarnings(glm(yv ~ gv, family = binomial()))
    co  <- summary(fit)$coef
    beta <- unname(co["gv","Estimate"])
    se   <- unname(co["gv","Std. Error"])
    zval <- unname(co["gv","z value"])
    pval <- unname(co["gv","Pr(>|z|)"])
    OR   <- exp(beta)
    c(zval, OR, pval)
  }, error = function(e) c(NA_real_, NA_real_, NA_real_))
  res
}

## ========= Parallel preparation =========
suppressWarnings({
  if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
})
library(doParallel)
library(data.table)

cl <- parallel::makeCluster(ncores)
registerDoParallel(cl)

## Divide the total number of times evenly among the cores
per_core <- as.integer(B_total %/% ncores)
rem      <- as.integer(B_total %% ncores)
task_vec <- rep(per_core, ncores)
if (rem > 0) task_vec[seq_len(rem)] <- task_vec[seq_len(rem)] + 1L

## ========= Parallel execution: one output file per core, written in batches =========
out_files <- foreach(core_id = seq_len(ncores), .packages = c("data.table")) %dopar% {
  Bi   <- task_vec[core_id]
  if (Bi == 0L) return(NA_character_)
  # Output file per core
  ofile <- file.path(out_dir,
                     sprintf("perm_core_%02d.csv%s", core_id, ifelse(use_gzip, ".gz", "")))
  # Write header (first time only)
  header_dt <- data.table(iter_global = integer(), Z = double(), OR = double(), P_add = double())
  if (use_gzip) {
    fwrite(header_dt, ofile, append = FALSE)
  } else {
    fwrite(header_dt, ofile, append = FALSE)
  }

  # The kernel's starting iteration number (optional)
  iter_start <- sum(task_vec[seq_len(core_id - 1L)], na.rm = TRUE)

  # Main loop: batch generation, permutation, fitting, and disk writing
  n_batches <- as.integer(ceiling(Bi / batch_size))
  idx_template <- seq_len(n)  # For shuffling
  wrote <- 0L

  for (b in seq_len(n_batches)) {
    this_batch <- min(batch_size, Bi - wrote)
    if (this_batch <= 0) break

    # Preallocation
    Z  <- numeric(this_batch)
    OR <- numeric(this_batch)
    P  <- numeric(this_batch)

    # In-batch loop
    for (k in seq_len(this_batch)) {
      # Shuffle y (permutate completely)
      yperm <- y[sample.int(n, n, replace = FALSE)]
      stats <- fit_once(yperm, g)
      Z[k]  <- stats[1]; OR[k] <- stats[2]; P[k] <- stats[3]
    }

    # Assemble and write to disk (append)
    iters <- (iter_start + wrote + seq_len(this_batch))
    dt <- data.table(iter_global = iters, Z = Z, OR = OR, P_add = P)
    fwrite(dt, ofile, append = TRUE)

    wrote <- wrote + this_batch
  }

  ofile
}

parallel::stopCluster(cl)

cat("Done. Per-core output files:\n")
print(out_files)
cat(sprintf("Total: %d replacements. The results are scattered in the directory: %s\n", 
            sum(task_vec), out_dir))





