#' re_training stochastic gradient boosting model with your data
#' args1 = new dataset (*.dummies.rmdup.pon.score)
#' args2 = output directory
#' args3 = prefix for output file
#' args4 = (optional) percent, 0-100%; the percent of the training set ( calls) will be used 
#' default is using 160 training samples

suppressMessages(library(tidyverse))
suppressMessages(library(gbm))
suppressMessages(library(ulimit))
ulimit::memory_limit(196000)

options(warn = -1)

args <- commandArgs(trailingOnly = F)
srcdir <- args[which(grepl('--file=', args))] %>% gsub('--file=', '', .) %>% dirname()
args <- args[(which(grepl('--args', args))+1):length(args)]

dir.create(file.path(args[2]), showWarnings = FALSE)

if (length(args) == 3){
  f <- read_tsv(paste0(dirname(srcdir), '/data/training_set_sampled_160.tsv.gz'), col_types = cols(.default = 'd', true_label = 'l')) %>%
    select(-c('delly', 'brass', 'dranger', 'snowman'))
  message('# SV calls from 160 PCAWG samples will be used together for re-training')

} else if (length(args) == 4){
  f <- read_tsv(paste0(dirname(srcdir), '/data/training_set_1711.tsv.gz'), col_types = cols(.default = 'd', true_label = 'l')) %>%
    select(-c('delly', 'brass', 'dranger', 'snowman'))
  set.seed(530)
  f <- f[sample(1:nrow(f), floor(nrow(f)*as.numeric(args[4])/100), replace = F), ]
  message(paste0('# ', args[4], '% of the training dataset will be used for re-training'))
  
} else {
  message('# ERROR: wrong inputs')
  quit()
}

# feature id
feature <- read_tsv(paste0(srcdir, '/feature_id.tsv'))
feature_col <- which(colnames(f) %in% feature$feature_prev)
colnames(f)[feature_col] <- feature$feature[order(factor(feature$feature_prev, levels = colnames(f)[feature_col]))]

# new dataset
g <- read_tsv(args[1], col_types = cols(.default = 'd', X1 = 'c', true_label = 'l')) %>%
  select(colnames(f))
m <- rbind(f, g)
rm(f); rm(g)

write_tsv(m, paste0(args[2], '/', args[3], '.training_set.tsv'))

# model training
set.seed(530)
m$true_label <- ifelse(m$true_label == T, 1, 0)
model.fit <- gbm(true_label ~ ., data = m, distribution = 'bernoulli', cv.folds = 5, n.trees = 10000, shrinkage = 0.01, interaction.depth = 6)

write_rds(model.fit, paste0(args[2], '/', args[3], '_gbm.fit.rds'))

