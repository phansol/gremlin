#' args1 = feature-annotated vcf (*.gremlin.feature.dummies.pon.score)
#' args2 = re-trained model.rds
#' args3 = (optional) threshold, 0-1

suppressMessages(library(tidyverse))
suppressMessages(library(gbm))
suppressMessages(library(data.table))

options(warn = -1)

args <- commandArgs(trailingOnly = T)

f <- read_tsv(args[1], col_types = cols(.default = 'd', X1 = 'c'))
f <- remove_rownames(f) %>% 
  column_to_rownames(var = 'X1')

model.fit <- readRDS(args[2])
best_iteration <- gbm.perf(model.fit, method = 'cv')
prob <- predict(model.fit, f, n.trees = best_iteration, type = 'response')
f$score_retrain <- round(prob, 5)

fwrite(x = f, file = paste0(args[1], '.retrained_score'), sep = '\t', row.names = T, col.names = T)

if (length(args) == 3){
  threshold <- as.numeric(args[3])
  
  f <- f %>% filter(score_retrain >= threshold)
  f$INFO <- rownames(f)
  
  f$SVTYPE <- ifelse(f$sv_type_DEL == 1, 'DEL', 
                     ifelse(f$sv_type_DUP == 1, 'DUP',
                            ifelse(f$sv_type_INV == 1, 'INV',
                                   ifelse(f$sv_type_TRA == 1, 'TRA','NA'))))
  
  f$tumor_depth_change_bpt1 <- round(f$tumor_depth_change_bpt1, 4) 
  f$tumor_depth_change_bpt2 <- round(f$tumor_depth_change_bpt2, 4) 
  f$tumor_var_mapq_bpt1 <- round(f$tumor_var_mapq_bpt1, 2) 
  f$tumor_var_mapq_bpt2 <- round(f$tumor_var_mapq_bpt2, 2)
  
  f <- f %>%
    separate('INFO', c('#CHR1', 'POS1', 'CHR2', 'POS2', 'CT'), sep = '_') %>%
    select(c('#CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE', 
             'tumor_var_read', 'tumor_var_split_read', 'tumor_var_sa_read',
             'normal_var_read', 'normal_var_split_read', 'normal_same_clip',
             'tumor_ref_read_bpt1', 'tumor_ref_read_bpt2', 'p_value_ref_var_read',
             'tumor_vaf_bpt1', 'tumor_vaf_bpt2', 
             'tumor_var_mapq_bpt1', 'tumor_var_mapq_bpt2',
             'tumor_depth_change_bpt1', 'tumor_depth_change_bpt2',
             'normal_other_var_cluster_bpt1', 'normal_other_var_cluster_bpt2',
             'normal_depth_bpt1', 'normal_depth_bpt2',
             'normal_panel', 'score', 'score_retrain'))
  
  write_tsv(f, paste0(args[1] %>% gsub('.gremlin.feature.dummies.pon.score', '', .), '.retrained_gremlin.somatic.svs.cutoff_', round(threshold, 3), '.tsv'))
}
