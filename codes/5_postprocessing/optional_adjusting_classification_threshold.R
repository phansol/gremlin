#' args1 = feature-annotated vcf (*.gremlin.dummies.rmdup.pon.score)
#' args2 = user-defined threshold (between 0-1)

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = T)

options(warn = -1)
f <- read_tsv(args[1], col_types = cols(.default = 'd', X1 = 'c')) %>%
  filter(score >= as.numeric(args[2]))

f$SVTYPE <- ifelse(f$sv_type_DEL == 1, 'DEL', 
                   ifelse(f$sv_type_DUP == 1, 'DUP',
                          ifelse(f$sv_type_INV == 1, 'INV',
                                 ifelse(f$sv_type_TRA == 1, 'TRA','NA'))))

f$tumor_depth_change_bpt1 <- round(f$tumor_depth_change_bpt1, 4) 
f$tumor_depth_change_bpt2 <- round(f$tumor_depth_change_bpt2, 4) 
f$tumor_var_mapq_bpt1 <- round(f$tumor_var_mapq_bpt1, 2) 
f$tumor_var_mapq_bpt2 <- round(f$tumor_var_mapq_bpt2, 2)

f <- f %>%
  separate('X1', c('#CHR1', 'POS1', 'CHR2', 'POS2', 'CT'), sep = '_') %>%
  select(c('#CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE', 
           'tumor_var_read', 'tumor_var_split_read', 'tumor_var_sa_read',
           'normal_var_read', 'normal_var_split_read', 'normal_same_clip',
           'tumor_ref_read_bpt1', 'tumor_ref_read_bpt2', 'p_value_ref_var_read',
           'tumor_vaf_bpt1', 'tumor_vaf_bpt2', 
           'tumor_var_mapq_bpt1', 'tumor_var_mapq_bpt2',
           'tumor_depth_change_bpt1', 'tumor_depth_change_bpt2',
           'normal_other_var_cluster_bpt1', 'normal_other_var_cluster_bpt2',
           'normal_depth_bpt1', 'normal_depth_bpt2',
           'normal_panel', 'score'))

write_tsv(f, args[1] %>% gsub('gremlin.feature.dummies.pon.score', paste0('gremlin.somatic.svs.cutoff_', round(as.numeric(args[2]), 3),'.tsv'), .))
