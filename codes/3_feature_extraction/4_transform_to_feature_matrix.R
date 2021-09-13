#' args[1] = vcf
#' args[2:5] = purity/ploidy/wgd_status/histology

suppressMessages(library(tidyverse))
suppressMessages(library(gsubfn))
options(warn = -1)

args <- commandArgs(trailingOnly = T)
options(scipen = 999)

log_file <- paste0(dirname(args[1]), '/log/', basename(args[1]) %>% gsub('.vcf.sort.SVinfo.fi.gc_mask', '', .), '.s4.success')
out_file <- args[1] %>% gsub('gremlin.sort.SVinfo.fi.gc_mask', 'gremlin.feature', .)

if (!(file.exists(log_file) & file.exists(out_file))){
  
  input <- read_tsv(args[1], col_types = cols(.default = 'c'))
  
  # parsing
  feature_matrix <- input[, c('MH', 'SVtype', 're_chr1', 're_pos1', 're_chr2', 're_pos2', 'terminal',
                              'T_ref1_count','T_ref2_count','T_all_discordant_frag','T_split_frag','T_sa_tag_frag', 'vaf_BP1', 'vaf_BP2',
                              'N_ref1_count', 'N_ref2_count', 'N_all_discordant_frag', 'N_split_frag', 'N_depth1', 'N_depth2',
                              'new_mate1', 'neo_mate1', 'new_mate2', 'neo_mate2', 'PairNormalSameClip',
                              'T_BP1_clip_readN', 'T_BP2_clip_readN', 'N_BP1_clip_readN', 'N_BP2_clip_readN',
                              'T_BP1_other_discordant_cluster', 'T_BP2_other_discordant_cluster',
                              'N_BP1_other_discordant_cluster', 'N_BP2_other_discordant_cluster',
                              'MAPQ1_min;med;max', 'MAPQ2_min;med;max', 'POS1_min;med;max', 'POS2_min;med;max',
                              'depth_ratio_change_bp1', 'depth_ratio_change_bp2',
                              'gc_content1', 'gc_content2', 'soft_masked1', 'soft_masked2')] 
  colnames(feature_matrix)[which(colnames(feature_matrix) == 'PairNormalSameClip')] <- 'N_same_clip'
  
  feature_matrix$vaf_BP1 <- gsub('%', '', feature_matrix$vaf_BP1) 
  feature_matrix$vaf_BP2 <- gsub('%', '', feature_matrix$vaf_BP2)
  
  # microhomology Y/N, the total number of new/neo mates at BP1 and BP2 --- if they exist, then the call may be true with a high probability
  col_to_numeric <- c('re_pos1', 're_pos2', 'vaf_BP1', 'vaf_BP2', 'T_ref1_count', 'T_ref2_count', 'T_all_discordant_frag', 'N_all_discordant_frag', 'N_ref1_count', 'N_ref2_count')
  if (nrow(feature_matrix) == 1){
    feature_matrix[col_to_numeric] <- sapply(feature_matrix[col_to_numeric], as.numeric) %>% as.list()
  } else {
    feature_matrix[col_to_numeric] <- sapply(feature_matrix[col_to_numeric], as.numeric)
  }
  
  for (i in 1:nrow(feature_matrix)){
    pos1 <- feature_matrix$re_pos1[i]
    pos2 <- feature_matrix$re_pos2[i]
    feature_matrix$distance[i] <- ifelse(feature_matrix$SVtype[i] == 'TRA', '250000000', as.integer(pos2 - pos1))  # distance b/w two break points, if SVtype=TRA(BND) then dist:=(the size of chr1)
    
    feature_matrix$MH[i] <- ifelse(feature_matrix$MH[i] == '.', F, T)    
    
    feature_matrix$new_mate1[i] <- gsubfn::strapply(feature_matrix$new_mate1[i], "\\(([^)]*)\\)", back = -1)[[1]] %>% as.numeric() %>% sum()
    feature_matrix$neo_mate1[i] <- gsubfn::strapply(feature_matrix$neo_mate1[i], "\\(([^)]*)\\)", back = -1)[[1]] %>% as.numeric() %>% sum()
    feature_matrix$new_mate2[i] <- gsubfn::strapply(feature_matrix$new_mate2[i], "\\(([^)]*)\\)", back = -1)[[1]] %>% as.numeric() %>% sum()
    feature_matrix$neo_mate2[i] <- gsubfn::strapply(feature_matrix$neo_mate2[i], "\\(([^)]*)\\)", back = -1)[[1]] %>% as.numeric() %>% sum()
    
    feature_matrix$vaf_BP1[i] <- ifelse(is.na(feature_matrix$vaf_BP1[i]), 0, feature_matrix$vaf_BP1[i])
    feature_matrix$vaf_BP2[i] <- ifelse(is.na(feature_matrix$vaf_BP2[i]), 0, feature_matrix$vaf_BP2[i])
    
    feature_matrix$`POS1_min;med;max`[i] <- ifelse(is.na(feature_matrix$`POS1_min;med;max`[i]), '0;1;2', feature_matrix$`POS1_min;med;max`[i])
    feature_matrix$`POS2_min;med;max`[i] <- ifelse(is.na(feature_matrix$`POS2_min;med;max`[i]), '0;1;2', feature_matrix$`POS2_min;med;max`[i])
    
    contingency_table <- rbind(c((feature_matrix$T_ref1_count[i] + feature_matrix$T_ref2_count[i]), 2*(feature_matrix$T_all_discordant_frag[i])), 
                               c((feature_matrix$N_ref1_count[i] + feature_matrix$N_ref2_count[i]), 2*(feature_matrix$N_all_discordant_frag[i]))) 
    feature_matrix$p_value[i] <- fisher.test(contingency_table, alternative = 'less')$p.value %>% round(., 6) %>% as.character() # p-value of Fisher's exact test
  }
  
  # the number of clusters of other discordant reads near the called break points --- the measure of noisiness at the break points
  feature_matrix$T_BP1_clip_readN <- as.integer(feature_matrix$T_BP1_clip_readN) - as.integer(feature_matrix$T_split_frag)
  feature_matrix$T_BP2_clip_readN <- as.integer(feature_matrix$T_BP2_clip_readN) - as.integer(feature_matrix$T_split_frag)
  
  feature_matrix$N_BP1_clip_readN <- as.integer(feature_matrix$N_BP1_clip_readN) - as.integer(feature_matrix$N_split_frag)
  feature_matrix$N_BP2_clip_readN <- as.integer(feature_matrix$N_BP2_clip_readN) - as.integer(feature_matrix$N_split_frag)
  
  # median MAPQ
  feature_matrix <- feature_matrix %>%
    separate(., `MAPQ1_min;med;max`, c('a','MAPQ1_med', 'b'), sep = ';') %>%
    dplyr::select(., -c('a','b'))
  feature_matrix <- feature_matrix %>%
    separate(., `MAPQ2_min;med;max`, c('a','MAPQ2_med', 'b'), sep = ';') %>%
    dplyr::select(., -c('a','b'))
  
  feature_matrix$MAPQ1_med[is.na(feature_matrix$MAPQ1_med)] <- mean(as.numeric(feature_matrix$MAPQ1_med), na.rm = T)
  feature_matrix$MAPQ2_med[is.na(feature_matrix$MAPQ2_med)] <- mean(as.numeric(feature_matrix$MAPQ2_med), na.rm = T)
  
  # the variation of reads' position T/F --- if the starting postion of all supporting reads are the same, it is highly likely false.
  feature_matrix$`POS1_min;med;max` <- feature_matrix$`POS1_min;med;max`=='0;0.0;0' 
  feature_matrix$`POS2_min;med;max` <- feature_matrix$`POS2_min;med;max`=='0;0.0;0' 
  colnames(feature_matrix)[which(colnames(feature_matrix)=='POS1_min;med;max')] <- 'POS1_variation'
  colnames(feature_matrix)[which(colnames(feature_matrix)=='POS2_min;med;max')] <- 'POS2_variation'
  
  # purity/ploidy/wgd_status/tumor_histology
  feature_matrix$purity <- args[2]
  feature_matrix$ploidy <- args[3]
  feature_matrix$wgd_status <- ifelse(args[4] == 'wgd', T, F)
  feature_matrix$histology <- args[5]
  
  # save
  colnames(feature_matrix)[colnames(feature_matrix) == 're_chr1'] <- 'CHR1'
  colnames(feature_matrix)[colnames(feature_matrix) == 're_chr2'] <- 'CHR2'
  colnames(feature_matrix)[colnames(feature_matrix) == 're_pos1'] <- 'POS1'
  colnames(feature_matrix)[colnames(feature_matrix) == 're_pos2'] <- 'POS2'
  colnames(feature_matrix)[colnames(feature_matrix) == 'terminal'] <- 'CT'
  
  write_tsv(feature_matrix, out_file)
  
}
