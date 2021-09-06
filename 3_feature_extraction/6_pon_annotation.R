#' sv panel of normal annotation
#' args1 = dummies
#' args2 = pon path

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

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


log_file <- paste0(dirname(args[1]), '/log/', basename(args[1]) %>% gsub('.feature.dummies', '', .), '.s6.success')
out_file <- paste0(args[1], '.pon')

if (!(file.exists(log_file) & file.exists(out_file))){
  
  pon_paths <- dir(args[2], pattern = 'panel_of_normals', full.names = T)
  pondel_path <- pon_paths[grepl('DEL', pon_paths)]
  pondup_path <- pon_paths[grepl('DUP', pon_paths)]
  poninv_path <- pon_paths[grepl('INV', pon_paths)]
  pontra_path <- pon_paths[grepl('TRA', pon_paths)]
  
  # input
  f <- read_tsv(args[1], col_types = cols(.default = 'd', X1 ='c', true_label = 'l')) %>%
    separate('X1', c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT'), sep = '_') 
  f$POS1 <- as.numeric(f$POS1)
  f$POS2 <- as.numeric(f$POS2)
  
  f$pon_sample <- 0
  fdel <- f %>% filter(SVtype_DEL == 1)
  fdup <- f %>% filter(SVtype_DUP == 1)
  finv <- f %>% filter(SVtype_INV == 1)
  ftra <- f %>% filter(SVtype_TRA == 1)
  
  # DEL
  if (nrow(fdel) != 0){
    pondel <- read_tsv(pondel_path, col_names = c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE', 'SAMPLE'), col_types = cols(POS1 = 'd', POS2 = 'd', .default = 'c'))
    pondel <- pon_match_prefix(fdel, pondel)
    for (i in 1:nrow(fdel)){
      sv_size <- fdel$POS2[i] - fdel$POS1[i]
      if (sv_size < 800){ # short sv
        pondel_i <- pondel %>% 
          filter(CHR1 == fdel$CHR1[i] & CHR2 == fdel$CHR2[i] & CT == fdel$CT[i] & abs(POS1 - fdel$POS1[i]) < sv_size/4 & abs(POS2 - fdel$POS2[i]) < sv_size/4)
      } else {
        pondel_i <- pondel %>% 
          filter(CHR1 == fdel$CHR1[i] & CHR2 == fdel$CHR2[i] & CT == fdel$CT[i] & abs(POS1 - fdel$POS1[i]) < dist & abs(POS2 - fdel$POS2[i]) < dist)
      }
      fdel$pon_sample[i] <- pondel_i$SAMPLE %>% unique() %>% length()
    }
    rm(pondel)
    write_tsv(fdel, paste0(args[1], '.del.pon.tmp'))
  }
  
  
  # DUP
  if (nrow(fdup) != 0){
    pondup <- read_tsv(pondup_path, col_names = c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE', 'SAMPLE'), col_types = cols(POS1 = 'd', POS2 = 'd', .default = 'c'))
    pondup <- pon_match_prefix(fdup, pondup)
    for (i in 1:nrow(fdup)){
      sv_size <- fdup$POS2[i] - fdup$POS1[i]
      if (sv_size < 800){ # short sv
        pondup_i <- pondup %>% 
          filter(CHR1 == fdup$CHR1[i] & CHR2 == fdup$CHR2[i] & CT == fdup$CT[i] & abs(POS1 - fdup$POS1[i]) < sv_size/4 & abs(POS2 - fdup$POS2[i]) < sv_size/4)
      } else {
        pondup_i <- pondup %>% 
          filter(CHR1 == fdup$CHR1[i] & CHR2 == fdup$CHR2[i] & CT == fdup$CT[i] & abs(POS1 - fdup$POS1[i]) < dist & abs(POS2 - fdup$POS2[i]) < dist)
      }
      fdup$pon_sample[i] <- pondup_i$SAMPLE %>% unique() %>% length()
    }
    rm(pondup)
    write_tsv(fdup, paste0(args[1], '.dup.pon.tmp'))
  }
  
  
  # INV
  if (nrow(finv) != 0){
    poninv <- read_tsv(poninv_path, col_names = c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE', 'SAMPLE'), col_types = cols(POS1 = 'd', POS2 = 'd', .default = 'c'))
    poninv <- pon_match_prefix(finv, poninv)
    for (i in 1:nrow(finv)){
      sv_size <- finv$POS2[i] - finv$POS1[i]
      if (sv_size < 800){ # short sv
        poninv_i <- poninv %>% 
          filter(CHR1 == finv$CHR1[i] & CHR2 == finv$CHR2[i] & CT == finv$CT[i] & abs(POS1 - finv$POS1[i]) < sv_size/4 & abs(POS2 - finv$POS2[i]) < sv_size/4)
      } else {
        poninv_i <- poninv %>% 
          filter(CHR1 == finv$CHR1[i] & CHR2 == finv$CHR2[i] & CT == finv$CT[i] & abs(POS1 - finv$POS1[i]) < dist & abs(POS2 - finv$POS2[i]) < dist)
      }
      finv$pon_sample[i] <- poninv_i$SAMPLE %>% unique() %>% length()
    }
    rm(poninv)
    write_tsv(finv, paste0(args[1], '.inv.pon.tmp'))
  }
  
  
  # TRA
  if (nrow(ftra) != 0){
    pontra <- read_tsv(pontra_path, col_names = c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE', 'SAMPLE'), col_types = cols(POS1 = 'd', POS2 = 'd', .default = 'c'))
    pontra <- pon_match_prefix(ftra, pontra)
    for (i in 1:nrow(ftra)){
      sv_size <- ftra$POS2[i] - ftra$POS1[i]
      pontra_i <- pontra %>% 
        filter(CHR1 == ftra$CHR1[i] & CHR2 == ftra$CHR2[i] & CT == ftra$CT[i] & abs(POS1 - ftra$POS1[i]) < dist & abs(POS2 - ftra$POS2[i]) < dist)
      ftra$pon_sample[i] <- pontra_i$SAMPLE %>% unique() %>% length()
    }
    rm(pontra)
    write_tsv(ftra, paste0(args[1], '.tra.pon.tmp'))
  }
  
  
  # save
  f_final <- rbind(fdel, fdup) %>%
    rbind(., finv) %>%
    rbind(., ftra)
  
  f_final$chr1_index <- f_final$CHR1 %>% gsub('chr', '', .) %>% gsub('X', 23, .) %>% gsub('Y', 24, .) %>% as.numeric()
  f_final$chr2_index <- f_final$CHR2 %>% gsub('chr', '', .) %>% gsub('X', 23, .) %>% gsub('Y', 24, .) %>% as.numeric()
  f_final <- f_final[order(f_final$CHR1, f_final$POS1, f_final$CHR2, f_final$POS2), ]
  
  f_final <- f_final %>%
    select(-c('chr1_index', 'chr2_index')) %>% 
    unite('INFO', c('CHR1', 'POS1', 'CHR2', 'POS2', 'CT'), sep = '_') %>%
    remove_rownames() %>%
    column_to_rownames(var = 'INFO')
  
  fwrite(x = f_final, file = out_file, sep = '\t', row.names = T, col.names = T)
  system(paste0('rm ', args[1], '.*.pon.tmp'))
  
}
