#' sv panel of normal annotation
#' args1 = feature-annotated vcf (*.dummies.rmdup.pon.score)
#' args2 = cohort-specific pon path/to/file (tab-delimited, without colnames, chr1 pos1 chr2 pos2 ct svtype sample)
#' args3 = cohort id

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

options(warn = -1)

args <- commandArgs(trailingOnly = T)

dist <- 200

pon_match_prefix <- function(vcf, pon){
  if (grepl('chr', vcf$CHR1[1]) & !(grepl('chr', pon$CHR1[1]))){
    pon$CHR1 <- paste0('chr', pon$CHR1)
    pon$CHR2 <- paste0('chr', pon$CHR2)
  } else if (!grepl('chr', vcf$CHR1[1]) & (grepl('chr', pon$CHR1[1]))){
    pon$CHR1 <- gsub('chr', '', pon$CHR1)
    pon$CHR2 <- gsub('chr', '', pon$CHR2)
  }
  return(pon)
}

# input
f <- read_tsv(args[1], col_types = cols(.default = 'd', X1 ='c', true_label = 'l')) %>%
  separate('X1', c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT'), sep = '_') 
f$POS1 <- as.numeric(f$POS1)
f$POS2 <- as.numeric(f$POS2)

f$normal_panel_cohort <- 0

if (nrow(f) != 0){
  pon <- read_tsv(args[2], col_names = c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE', 'SAMPLE'), col_types = cols(POS1 = 'd', POS2 = 'd', .default = 'c'))
  pon <- pon_match_prefix(f, pon)
  for (i in 1:nrow(f)){
    sv_size <- f$POS2[i] - f$POS1[i]
    if (sv_size < 800 & SVtype_TRA == 0){ # short sv
      pon_i <- pon %>% 
        filter(CHR1 == f$CHR1[i] & CHR2 == f$CHR2[i] & CT == f$CT[i] & abs(POS1 - f$POS1[i]) < sv_size/4 & abs(POS2 - f$POS2[i]) < sv_size/4)
    } else {
      pon_i <- pon %>% 
        filter(CHR1 == f$CHR1[i] & CHR2 == f$CHR2[i] & CT == f$CT[i] & abs(POS1 - f$POS1[i]) < dist & abs(POS2 - f$POS2[i]) < dist)
    }
    f$normal_panel_cohort[i] <- pon_i$SAMPLE %>% unique() %>% length()
  }
}
colnames(f)[which(colnames(f) == 'normal_panel_cohort')] <- paste0('normal_panel_', args[3])

f <- f %>%
  unite('INFO', c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT'), sep = '_') %>%
  remove_rownames() %>%
  column_to_rownames(var = 'INFO')

fwrite(x = f, file = paste0(args[1], '.pon_', args[3]), sep = '\t', row.names = T, col.names = T)
