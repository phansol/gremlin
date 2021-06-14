#' sv filtering using gremlin
#' args1 = feature_matrix.dummies.pon

suppressMessages(library(tidyverse))
suppressMessages(library(gbm))
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = F)
srcdir <- args[which(grepl('--file=', args))] %>% gsub('--file=', '', .) %>% dirname()
args <- args[(which(grepl('--args', args))+1):length(args)]

outdir <- args[1] %>% dirname() %>% dirname()
sampleid <- args[1] %>% basename() %>% gsub('.feature.dummies.pon', '', .)

# GREMLIN
gremlin <- readRDS(paste0(srcdir, '/gremlin_gbm.fit.rds'))
best_iteration <- 10000

threshold_1 <- 0.885457596188905 # 95% sensitivity, 77% precision on the PCAWG validation set
threshold_2 <- 0.570132241490398 # 97.5% sensitivity, 51% precision on the PCAWG validation set

# input
f <- read_tsv(args[1], col_types = cols(.default = 'd', X1 = 'c'))
f <- remove_rownames(f) %>% 
  column_to_rownames(var = 'X1')

# classification
prob <- predict(gremlin, f, n.trees = best_iteration, type = 'response')
f$score <- round(prob, 5)

# feature id
feature <- read_tsv(paste0(srcdir, '/feature_id.tsv'))
feature_col <- which(colnames(f) %in% feature$feature_prev)
colnames(f)[feature_col] <- feature$feature[order(factor(feature$feature_prev, levels = colnames(f)[feature_col]))]
fwrite(x = f, file = paste0(outdir, '/', basename(args[1]), '.score'), sep = '\t', row.names = T, col.names = T)

# formatting
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
           'normal_panel', 'score'))

# save
tier1 <- f %>% filter(score >= threshold_1)
tier2 <- f %>% filter(score < threshold_1 & score >= threshold_2)

write_tsv(tier1, paste0(outdir, '/', sampleid, '.sv.gremlin.tier1.vcf'))
write_tsv(tier2, paste0(outdir, '/', sampleid, '.sv.gremlin.tier2.vcf'))
