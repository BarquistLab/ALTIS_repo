create_dfit <- function(data_table, sample_stats, lab, description, res_dir, lab_in="INPUT", lab_out="OUTPUT"){
  
  dfit <- pivot_longer(data_table, !c("seq", "pos", "locus_tag", "k", "type", "TI"), names_to = "sample", values_to = "raw")
  
  dfit_NZ <- subset(dfit, grepl(lab_in, sample) & raw>0)
  data_table <- subset(data_table, TI %in% dfit_NZ$TI)
  dfit <- subset(dfit, TI %in% dfit_NZ$TI)
  
  data_table <- data_table %>% arrange(grepl("neg", TI), seq, pos)
  data_table$locus_tag <- factor(data_table$locus_tag, levels=data_table$locus_tag %>% unique())
  data_table$it_g <- as.integer(data_table$locus_tag)
  read_counts <- data_table[, grep(paste(lab_in, lab_out, sep="|"), colnames(data_table))]
  dfit$locus_tag <- factor(dfit$locus_tag, levels=data_table$locus_tag %>% unique())
  dfit$it_g <- as.integer(dfit$locus_tag)
  
  # normalization with pseudo genes
  lib_sizes <- colSums(data_table[, grep(paste(lab_in, lab_out, sep="|"), colnames(data_table))])
  
  dfit$com_id <- data_table$com_id[match(dfit$locus_tag, data_table$locus_tag)]
  dfit$N <- lib_sizes[match(dfit$sample, names(lib_sizes))]
  
  dfit$cpm <- (dfit$raw+0.5)/(dfit$N+1)*1e+06
  dfit$logcpm <- log(dfit$cpm)
  
  dfit$exp_id <- sample_stats$exp_id[match(dfit$sample, sample_stats$sample)]
  dfit$TI_exp_id <- paste(dfit$TI, dfit$exp_id, sep=":")
  
  # calculate number of sites per gene
  non_zero_TI <- subset(dfit, grepl(lab_in, sample) & raw!=0)$TI_exp_id

  insertions <- data_table %>% group_by(it_g, locus_tag) %>% summarise(n_g=length(locus_tag))
  n_gti <- insertions$n_g
  N_ti <- sum(n_gti)
  
  # keep only IS non-zero input
  dfit_NZ <- subset(dfit, grepl(lab_in, sample) & raw!=0)
  dfit_NZ_out <- subset(dfit, grepl(lab_out, sample) & TI_exp_id %in% dfit_NZ$TI_exp_id) %>% dplyr::select(c("TI_exp_id", "raw", "N"))
  colnames(dfit_NZ_out) <- c("TI_exp_id", "raw_out", "N_out")
  dfit_NZ <- cbind(dfit_NZ, dfit_NZ_out[match(dfit_NZ$TI_exp_id, dfit_NZ_out$TI_exp_id),] %>% dplyr::select(c("raw_out", "N_out")))
  
  dfit_NZ <- dfit_NZ %>% mutate(alpha_in=log(raw))
  
  sample_stats2 <- dfit_NZ %>% group_by(exp_id, sample) %>% summarise(alpha_sa=mean(alpha_in))
  dfit_NZ$alpha_sa <- sample_stats2$alpha_sa[match(dfit_NZ$exp_id, sample_stats2$exp_id)]
  
  # save count table, processed (long) data frame (dfit_NZ) and append description to summary file
  fwrite(data_table, paste0(res_dir, "counts_", lab, ".csv"))
  fwrite(dfit_NZ, paste0(res_dir, "dfit_", lab, ".csv"))
  
  rm(list=ls()[grep("dfit", ls())])
}

# export TraDIS data, input-output pairs
export_cmdstan <- function(dfit_file, counts_file, cmd_dir, lab, lab_in="INPUT", lab_out="OUTPUT",
                           rm_genes=NULL, norm_fac, a_in=c(0,-0.7), b=0.2, organism="Salmonella", count_min=0){
  data_table <- fread(counts_file, data.table=F)
  
  rm_IS <- which(rowSums(data_table[, grep(paste(lab_in, lab_out, sep="|"), colnames(data_table))])<count_min)
  if(length(rm_IS)>0) data_table <- data_table[-rm_IS, ]
  if(!is.null(rm_genes)) data_table <- data_table[-which(data_table$locus_tag %in% rm_genes),]
  
  data_table$locus_tag <- factor(data_table$locus_tag, levels=data_table$locus_tag %>% unique())
  data_table$it_g <- as.integer(data_table$locus_tag)
  
  dfit_NZ <- fread(dfit_file, data.table=F)
  
  sample_stats <- subset(dfit_NZ, type=="pseudo") %>% group_by(exp_id, sample) %>%
    summarise(N_nz=length(which(raw_out>0)), n_g=length(raw),
              N_in=mean(N), N_out=mean(N_out))
  sample_stats$Z <- (sample_stats$n_g-sample_stats$N_nz)/sample_stats$n_g
  
  raw_counts_in <- t(data_table[, sample_stats$sample])
  raw_counts_out <- t(data_table[, sub(lab_in, lab_out, sample_stats$sample)])
  
  it_g <- data_table$it_g
  N_ti <- length(it_g)
  N_g <- max(it_g)
  
  Z_ps <- sample_stats$Z
  N_exp <- length(Z_ps)
  
  logFC_g <- rep(-10, N_g)
  
  if(!is.null(rm_genes)) dfit_NZ <- dfit_NZ %>% subset(!(locus_tag %in% rm_genes))
  gene_stats <- dfit_NZ %>% group_by(locus_tag, it_g) %>% summarise(n_gz=length(which(raw_out==0)), n_g=length(raw), logcpm_in=mean(logcpm))
  gene_stats$drop_out <- gene_stats$n_gz/gene_stats$n_g
  
  sample_stats <- dfit_NZ %>% group_by(exp_id, sample) %>% summarise(N_nz=length(which(raw_out>0)), n_g=length(raw), N_in=mean(N), N_out=mean(N_out))
  N_ti_max <- max(sample_stats$n_g)
  
  data_file <- paste(cmd_dir, organism, "_", lab, "_data.R", sep="")
  
  dfit_NZ$n_ps <- exp(norm_fac)[match(dfit_NZ$sample, sample_stats$sample)]
  dfit_NZ$alpha_in <- log(dfit_NZ$raw * dfit_NZ$n_ps)
  sample_stats2 <- dfit_NZ %>% group_by(exp_id, sample) %>% summarise(alpha_sa=mean(alpha_in))
  dfit_NZ$alpha_sa <- sample_stats2$alpha_sa[match(dfit_NZ$exp_id, sample_stats2$exp_id)]
  
  stan_rdump(c("N_ti", "N_g", "N_exp", "Z_ps", "it_g", "raw_counts_in", "raw_counts_out", "N_ti_max", "norm_fac", "b", "a_in"),file=data_file)
  
  init_file <- paste(cmd_dir, organism, "_", lab, "_init.R", sep="")
  stan_rdump(c("logFC_g"),file=init_file)
  
  rm(raw_counts_in, raw_counts_out, data_table)
  dfit_NZ
}

# extract sample information from sample name
sample_stats_shigella <- function(sample, it_sa){
  sample_df <- data.frame(sample, it_sa)
  sample_df$it_c <- strsplit(sample %>% as.character(), "_") %>% sapply(function(x) x[2]) %>% as.factor() %>% as.numeric()
  sample_df$exp_id <- strsplit(sample %>% as.character(), "_") %>% sapply(function(x) sub("Lib", "", x[1])) %>% as.numeric()
  sample_df
}
