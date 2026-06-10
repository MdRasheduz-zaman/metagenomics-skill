#!/usr/bin/env Rscript
# DADA2 ASV inference for one amplicon sample. Driven by metagx (dada2.yaml registry):
# user flags (--trunc_len_f/--trunc_len_r/--max_ee_f/--max_ee_r/--trunc_q/--pool) plus
# workflow-injected I/O (--r1/--r2/--out_table/--out_seqs/--threads). Paired if --r2 given.
#
# Outputs:
#   <out_table>: TSV with columns asv_id, <sample>, sequence (per-ASV read counts)
#   <out_seqs> : FASTA of the ASV sequences
suppressMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}

r1        <- getarg("--r1")
r2        <- getarg("--r2", "")
out_table <- getarg("--out_table")
out_seqs  <- getarg("--out_seqs")
sample    <- getarg("--sample", "sample")
threads   <- as.integer(getarg("--threads", "1"))
truncF    <- as.integer(getarg("--trunc_len_f", "0"))
truncR    <- as.integer(getarg("--trunc_len_r", "0"))
maxeeF    <- as.numeric(getarg("--max_ee_f", "2"))
maxeeR    <- as.numeric(getarg("--max_ee_r", "2"))
truncQ    <- as.integer(getarg("--trunc_q", "2"))
pool      <- getarg("--pool", "independent")
pool_val  <- if (pool == "true") TRUE else if (pool == "pseudo") "pseudo" else FALSE

paired <- nzchar(r2)
filtdir <- file.path(dirname(out_table), "dada2_filt")
dir.create(filtdir, showWarnings = FALSE, recursive = TRUE)
fF <- file.path(filtdir, "F_filt.fastq.gz")

if (paired) {
  fR <- file.path(filtdir, "R_filt.fastq.gz")
  filterAndTrim(r1, fF, r2, fR, truncLen = c(truncF, truncR),
                maxEE = c(maxeeF, maxeeR), truncQ = truncQ, rm.phix = TRUE,
                compress = TRUE, multithread = threads)
  errF <- learnErrors(fF, multithread = threads)
  errR <- learnErrors(fR, multithread = threads)
  ddF <- dada(fF, err = errF, pool = pool_val, multithread = threads)
  ddR <- dada(fR, err = errR, pool = pool_val, multithread = threads)
  merged <- mergePairs(ddF, fF, ddR, fR)
  seqtab <- makeSequenceTable(merged)
} else {
  filterAndTrim(r1, fF, truncLen = truncF, maxEE = maxeeF, truncQ = truncQ,
                rm.phix = TRUE, compress = TRUE, multithread = threads)
  errF <- learnErrors(fF, multithread = threads)
  ddF <- dada(fF, err = errF, pool = pool_val, multithread = threads)
  seqtab <- makeSequenceTable(ddF)
}

seqtab <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads)
seqs <- colnames(seqtab)
counts <- as.integer(seqtab[1, ])
ids <- sprintf("ASV%04d", seq_along(seqs))

writeLines(paste0(">", ids, "\n", seqs), out_seqs)
df <- data.frame(asv_id = ids, count = counts, sequence = seqs)
colnames(df)[2] <- sample
write.table(df, out_table, sep = "\t", quote = FALSE, row.names = FALSE)
