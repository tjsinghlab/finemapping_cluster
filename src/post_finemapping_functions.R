
##reads in results
read_in_finemap_res <- function(file_path, pattern){
  file_list = list.files(paste0("/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/", "/",file_path),
                         pattern = pattern, full.names = TRUE)
  print(length(file_list))
  dfM.finemapped.0 = lapply(file_list, fread)
  dfM.finemapped.0 = rbindlist(dfM.finemapped.0, fill = T)
  return(dfM.finemapped.0)
}

get_rsids_by_PIP_2.0 <- function(sumstats_name, ld_pop, window_mb, LDpanel){
  susie <- read_in_finemap_res(
    paste0("output/", sumstats_name, "/", ld_pop,
           "_", window_mb, "Mb_window/susie"), glob2rx(paste0(LDpanel, "*adjustedvar.tsv")))
  
  carma <- read_in_finemap_res(
    paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb_window/CARMA"), 
    glob2rx(paste0("*", LDpanel, "*.txt.gz")))
  
  susie <- susie %>% 
    rename_with(~ paste0("susie_", .x),starts_with("CS")) %>% 
    select(rsid, locus, prob = PIP, beta, se, A1 = allele1, A2 = allele2, chr = chromosome, bp = position, CSsusie= CS) %>% 
    mutate(method = "susie") 
  
  carma <- carma %>% 
    select(rsid, locus, prob = PIP, beta, se, A1 = allele1, A2 = allele2, chr = chromosome, bp = position, outlier, CScarma = CS) %>% 
    mutate(method = "carma")
  
  res <- rbind(susie, carma, fill = T)
  res$LDpanel <- LDpanel
  res$window <- window_mb
  return(res)
}