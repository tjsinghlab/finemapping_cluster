
generate_per_locus_data_and_submit_pipeline <- function(df.leadSNP, path_df_sumstats, sumstats_name, N_tot, N_cases, window_mb,  ld_pop){
  window_bp <- window_mb*1000000
  df.leadSNP <- df.leadSNP %>%
    mutate(LOWER = ifelse((BP-window_bp) < 0, 1, BP-window_bp), #if lead snp is too close to chr beginning
           UPPER = BP+window_bp)
  
  dir.create(paste0("output/", sumstats_name))
  dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb_window"))
  dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb_window/ss"))
  
  for (i in 1:nrow(df.leadSNP)) {
    LOCUS <- df.leadSNP$locus[i]
    CHR <- df.leadSNP$CHR[i]
    START <- df.leadSNP$LOWER[i]
    END <- df.leadSNP$UPPER[i]
    job_name <- paste("run_FM_pipeline", sumstats_name, LOCUS, paste0(window_mb, "Mb"), sep = "_")
    system(paste("sbatch -J", job_name,
                 "src/run_finemapping_per_locus.sh", 
                 sumstats_name, 
                 path_df_sumstats, 
                 ld_pop,
                 LOCUS,
                 CHR,
                 START,
                 END,
                 N_tot,
                 N_cases,
                 window_mb, sep = " "))
  }
}

run_CARMA <- function(sumstat, ld, sumstats_name, ld_pop, window_mb, locus, LDpanel){
  z.list<-list()
  ld.list<-list()
  lambda.list<-list()
  sumstat$Z <- sumstat$beta/sumstat$se
  z.list[[1]]<-sumstat$Z
  ld.list[[1]]<-as.matrix(ld)
  lambda.list[[1]]<-1 
  CARMA.results<-CARMA(z.list, ld.list, lambda.list=lambda.list, outlier.switch=TRUE, rho.index = 0.95)
  saveRDS(CARMA.results, file = paste0("output/",sumstats_name,"/", ld_pop, "_", 
                                       window_mb, "Mb_window/CARMA/", sumstats_name, "_", LDpanel, "_", 
                                       window_mb, "Mb_locus_",locus,".rds"))
  print("Outliers")
  print(CARMA.results$Outliers)
  
  ###### Posterior inclusion probability (PIP) and credible set (CS)
  sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0) 
  
  if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
    for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){ 
      sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
    } }
  
  sumstat.result$outlier <- 0
  sumstat.result$outlier[CARMA.results[[1]]$Outliers$Index] <- 1
  
  dir.create(paste0("output/",sumstats_name,"/", ld_pop, "_", 
                    window_mb, "Mb_window/CARMA/"))
  
  ###### write the GWAS summary statistics with PIP and CS
  fwrite(x = sumstat.result,
         file = paste0("output/",sumstats_name,"/", ld_pop, "_", 
                       window_mb, "Mb_window/CARMA/", sumstats_name, "_", LDpanel, "_", 
                       window_mb, "Mb_locus_",locus,".txt.gz"), 
         sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
  
}


run_susie <- function(sumstats, LDmat, N_tot, N_cases, sumstats_name, ld_pop,  window_mb, locus, LDpanel, coverage = 0.95){
  window_bp <- window_mb*1000000
  if(length(sumstats$rsid) != nrow(LDmat)){
    print("Problem! SS dimensions do not agree with LD matrix dimensions")
  }
  colnames(LDmat) <- rownames(LDmat) <- paste(sumstats$rsid,sumstats$allele1, sumstats$allele2, sep = ":")
  
  dir.create(paste0("output/",sumstats_name,"/",ld_pop,"_", 
                    window_mb, "Mb_window/susie_res_", window_mb, "Mb"))
  phi <- N_cases/N_tot
  
  susie_output <- susie_rss(R = as.matrix(LDmat), n = N_tot, vary_y = 1/(phi*(1-phi)), bhat = sumstats$beta, shat = sumstats$se)
  
  if(susie_output$converged == TRUE){
    print("converged!")
    sumstats$converged = susie_output$converged
    sumstats$PIP = susie_output$pip
    cs_out <- susie_output$sets
    if(!is.null(cs_out$cs)){
      print(locus)
      print(cs_out$cs)
      sumstats$num_credsets <- length(cs_out$cs)
      for (j in 1:length(cs_out$cs)) {
        #sumstats[cs_out$cs[[j]], "credibleSets"]  = paste0("L",j)
        sumstats[,paste0("CS", j)] <- 0
        sumstats[cs_out$cs[[j]], paste0("CS", j)] <- 1
      }
    }else{
      print("NUll!")
    }
    #susie_plot(susie_output, y = "PIP", add_bar = TRUE, add_legend = TRUE, main = LOCUS)
    
    ## write finemapped sumstats file
    write_tsv(sumstats, paste0("output/",sumstats_name,"/",ld_pop,"_", 
                               window_mb, "Mb_window/susie_res_", 
                               window_mb, "Mb/", LDpanel,"_locus", locus, "_cov", coverage,"_adjustedvar.tsv"))
  }
}


get_ss_per_locus <- function(path_df_sumstats, CHR, LOCUS, START, END){
  df.sumstats <- fread(path_df_sumstats)
  df_chrom <- df.sumstats %>% 
    filter(chromosome == CHR) %>% 
    arrange(position)
  ROW_START = which(df_chrom[,"position"]>START)[1]
  ROW_END = which(df_chrom[,"position"]<END)[length(which(df_chrom[,"position"]<END))]
  df_chrom[ROW_START:ROW_END,"locus"] = LOCUS
  return(df_chrom[ROW_START:ROW_END,])
}


get_ld_per_locus <- function(ss_per_locus, LOCUS, CHR, START, END){
  path_to_LD <- find_LD_block(CHR, START, END)
  print(path_to_LD)
  bim <- fread(paste0(path_to_LD, ".bim"))
  ld <- fread(paste0(path_to_LD, ".ld"))
  colnames(bim) <- c("chr", "rsid", "dk", "pos", "alt", "ref")
  ss_per_locus$SNP1 <- paste(ss_per_locus$chromosome,  #allele1 = alt, no need to change direction of beta
                             ss_per_locus$position, 
                             ss_per_locus$allele1, 
                             ss_per_locus$allele2, sep = ":")
  ss_per_locus$SNP2 <- paste(ss_per_locus$chromosome, 
                             ss_per_locus$position, 
                             ss_per_locus$allele2, 
                             ss_per_locus$allele1, sep = ":")
  bim$SNP <- paste(bim$chr, bim$pos, bim$alt, bim$ref, sep = ":")
  bim$ss_matched <- bim$SNP %in% c(ss_per_locus$SNP2, ss_per_locus$SNP1)
  ld_filtered <- as.matrix(ld)[bim$ss_matched, bim$ss_matched]
  if(sum(is.na(ld_filtered)) != 0){
    print("LD mat contains NAs")
  }else{
    bim_filtered <- bim %>% filter(ss_matched)
    ss_filtered <- rbind(ss_per_locus %>% 
                           filter(SNP1 %in% bim_filtered$SNP) %>% 
                           rename(SNP = SNP1) %>% 
                           select(!c(SNP2)),
                         ss_per_locus %>% 
                           filter(SNP2 %in% bim_filtered$SNP) %>% 
                           mutate(beta = -1*beta)  %>% #converts beta to be in same direction as ld allele
                           rename(SNP = SNP2, allele1 = allele2, allele2 = allele1) %>%  #then we rename alleles
                           select(!c(SNP1)))
    
    ss_filtered <- inner_join(bim %>% select(SNP), ss_filtered, by = "SNP")
    print(paste(sum(ss_filtered$SNP != bim_filtered$SNP), "variants out of order")) #check vars are in correct order
    return(list(ss_filtered, ld_filtered))
  }
}


find_LD_block <- function(CHR, START, END){
  approx_blocks <- fread("/gpfs/commons/groups/nygcfaculty/sghatan/UKb_LDmatrices/approx_LD_blocks.txt")
  if(nrow(approx_blocks %>% filter(chr == CHR, start < START, stop > END)) == 1){
    LD_start <- approx_blocks %>% filter(chr == CHR, start < START, stop > END) %>% pull(start)
    LD_stop <- approx_blocks %>% filter(chr == CHR, start < START, stop > END) %>% pull(stop)
    filepath <- paste0("/gpfs/commons/groups/nygcfaculty/sghatan/UKb_LDmatrices/chr", CHR, "/", LD_start, ".", LD_stop, "/", LD_start, ".", LD_stop)
    return(filepath)
  }else{
    print("locus overlaps multiple LD blocks!") #make exception
  }
}

