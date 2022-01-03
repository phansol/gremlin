#' remove duplicates & pre-filtering
#' calls with re-annotated T_all_discordant_frag <= 1 OR N_sa_tag_frag > 0 are false.
#' args1 = vcf
#' args2 = sampleid 

args <- commandArgs(trailingOnly = T)
options(warn = -1)

suppressMessages(library(tidyverse))

outdir <- dirname(args[1])
sampleid <- args[2]

log_file <- paste0(outdir, '/log/', sampleid, '.s2.success')
out_file <- paste0(args[1], '.fi')

if (!(file.exists(log_file) & file.exists(out_file))){
  
  vcf <- read_tsv(args[1], col_types = cols(.default = 'c'))
  
  # remove duplicates
  vcf_reassigned <- vcf %>% select(., c('re_chr1', 're_pos1', 're_chr2', 're_pos2', 'terminal'))
  vcf_rmdup <- vcf[!duplicated(vcf_reassigned), ] 
  
  # prefiltering
  vcf_rmdup <- vcf_rmdup %>% 
    tidyr::separate(., `Tumor_Ref1;Ref2;AllDiscordantFragments;SplitFragments;SATagFragments;Vaf1;Vaf2`, 
                    c('T_ref1_count','T_ref2_count','T_all_discordant_frag','T_split_frag','T_sa_tag_frag', 'vaf_BP1', 'vaf_BP2'), sep = ';') %>% 
    tidyr::separate(., `PairNormal_Ref1;Ref2;AllDiscordantFragments;SplitFragments;SATagFragments;FragCount1;FragCount2`, 
                    c('N_ref1_count', 'N_ref2_count', 'N_all_discordant_frag', 'N_split_frag', 'N_sa_tag_frag', 'N_depth1', 'N_depth2'), sep = ';') 
  
  vcf_rmdup$T_all_discordant_frag <- as.numeric(vcf_rmdup$T_all_discordant_frag)
  vcf_rmdup$N_sa_tag_frag <- as.numeric(vcf_rmdup$N_sa_tag_frag)
  
  vcf_filtered <- vcf_rmdup %>% dplyr::filter(T_all_discordant_frag >= 2 & N_sa_tag_frag == 0) 
  
  if (nrow(vcf_filtered) == 0){
    writeLines(c('# Prefiltered vcf is empty', '# Job finished'), paste0(outdir, '/log/', sampleid, '.s2.job.finished'))
  } else {
    write_tsv(vcf_filtered, out_file)
  }
  
}
