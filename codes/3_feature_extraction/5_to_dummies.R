#' args1 = *.gremlin.feature

suppressMessages(library(tidyverse))
suppressMessages(library(dummies))
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = T)

log_file <- paste0(dirname(args[1]), '/log/', basename(args[1]) %>% gsub('.feature', '', .), '.s5.success')
out_file <- paste0(args[1], '.dummies')

if (!(file.exists(log_file) & file.exists(out_file))){
  
  f <- read_tsv(args[1], col_types = cols(.default = 'c')) %>%
    tidyr::unite(., 'INFO', c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT'), sep = '_') 
  
  f <- remove_rownames(f) %>% 
    column_to_rownames(var = 'INFO')
  
  # Transform categorical variables to dummy variables
  f <- dummy.data.frame(f, names = c('MH', 'SVtype', 'POS1_variation', 'POS2_variation', 'wgd_status', 'histology'), sep= "_", omit.constants = F)
  dummy_cols <- c(MH_FALSE = 0, MH_TRUE = 0, 
                  SVtype_DEL = 0, SVtype_DUP = 0, SVtype_INV = 0, SVtype_TRA = 0, 
                  POS1_variation_TRUE = 0, POS1_variation_FALSE = 0, 
                  POS2_variation_TRUE = 0, POS2_variation_FALSE = 0, 
                  wgd_status_TRUE = 0, wgd_status_FALSE = 0,
                  histology_Biliary = 0, histology_Bladder = 0, histology_Bone_SoftTissue = 0, histology_Breast = 0, 
                  histology_Cervix = 0, histology_CNS = 0, histology_Colon_Rectum = 0, histology_Esophagus = 0, 
                  histology_Head_Neck = 0, histology_Hematologic = 0, histology_Kideny = 0, histology_Liver = 0, 
                  histology_Lung = 0, histology_Ovary = 0, histology_Pancreas = 0, histology_Prostate = 0, 
                  histology_Skin = 0, histology_Stomach = 0, histology_Thyroid = 0, histology_Uterus = 0)
  f <- add_column(f, !!!dummy_cols[setdiff(names(dummy_cols), names(f))]) 
  
  col_order <- c('T_ref1_count', 'T_ref2_count', 'T_all_discordant_frag', 'T_split_frag', 'T_sa_tag_frag', 'vaf_BP1', 'vaf_BP2',
                 'N_ref1_count', 'N_ref2_count', 'N_all_discordant_frag', 'N_split_frag', 'N_depth1', 'N_depth2', 
                 'new_mate1', 'neo_mate1', 'new_mate2', 'neo_mate2', 'N_same_clip', 
                 'T_BP1_clip_readN', 'T_BP2_clip_readN', 'N_BP1_clip_readN', 'N_BP2_clip_readN', 
                 'T_BP1_other_discordant_cluster', 'T_BP2_other_discordant_cluster', 
                 'N_BP1_other_discordant_cluster', 'N_BP2_other_discordant_cluster', 
                 'MAPQ1_med', 'MAPQ2_med', 'depth_ratio_change_bp1', 'depth_ratio_change_bp2', 
                 'gc_content1', 'gc_content2', 'soft_masked1', 'soft_masked2', 
                 'distance', 'p_value', 'purity', 'ploidy', 
                 'MH_FALSE', 'MH_TRUE',   'SVtype_DEL', 'SVtype_DUP', 'SVtype_INV', 'SVtype_TRA', 
                 'POS1_variation_TRUE', 'POS1_variation_FALSE', 'POS2_variation_TRUE', 'POS2_variation_FALSE', 
                 'wgd_status_TRUE', 'wgd_status_FALSE', 
                 'histology_Biliary', 'histology_Bladder', 'histology_Bone_SoftTissue', 'histology_Breast', 
                 'histology_Cervix', 'histology_CNS', 'histology_Colon_Rectum', 'histology_Esophagus', 
                 'histology_Head_Neck', 'histology_Hematologic', 'histology_Kideny', 'histology_Liver', 
                 'histology_Lung', 'histology_Ovary', 'histology_Pancreas', 'histology_Prostate', 
                 'histology_Skin', 'histology_Stomach', 'histology_Thyroid', 'histology_Uterus')
  
  
  f <- f %>% dplyr::select(col_order)
  
  fwrite(x = f, file = out_file, sep = '\t', col.names = T, row.names = T)
  
}
