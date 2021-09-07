check_unknown <- function(string){
  d <- grepl('decoy', string)
  r <- grepl('random', string)
  u <- grepl('chrUn', string)
  gl <- grepl('GL', string)
  nc <- grepl('NC', string)
  hs <- grepl('hs37', string)
  return(d + r + u + gl + nc + hs)
}

chr_to_int <- function(x){
  y <- x %>% gsub('chr', '', .) %>% gsub('X', 23, .) %>% gsub('Y', 24, .) %>% as.numeric()
  return(y)
}

split_pos <- function(x){
  return(substr(x, 1, nchar(x) - 1))
}

split_dir <- function(x){
  dir <- substr(x, nchar(x), nchar(x))
  dir <- ifelse(dir == '[', 5, 3)
  return(dir)
}

vcf_formatting <- function(vcf_path, ref.fai_path){
  primary_chr <- c(1:22, 'X', 'Y', paste0('chr', c(1:22, 'X', 'Y')))
  
  # input
  vcf <- read_tsv(vcf_path, col_names = c('CHR1', 'POS1', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR'), comment = '#',
                  col_types = cols('CHR1' = 'c')) %>% 
    separate(ALT, into = c('chr2', 'pos2'), sep = ':') %>%
    separate(chr2, into = c('a', 'CHR2'), sep = '\\[|\\]') %>%
    separate(INFO, into = c('b', 'c'), sep = 'MATEID=', remove = F) %>%
    separate(c, into = c('mate', 'd'), sep = ';') %>%
    select(c('CHR1', 'POS1', 'CHR2', 'pos2', 'ID', 'mate')) %>%
    filter(CHR1 %in% primary_chr & CHR2 %in% primary_chr)
  
  if (nrow(vcf) == 0){
    return('empty')
  } else {
    vcf$pos2 <- gsub('([[:alpha:]])', '', vcf$pos2)
    vcf$POS2 <- sapply(vcf$pos2, split_pos)
    vcf$dir2 <- sapply(vcf$pos2, split_dir)
    
    for (i in 1:nrow(vcf)){
      mate <- vcf %>% filter(ID == vcf$mate[i])
      vcf$dir1[i] <- mate$dir2
    }
    
    vcf <- vcf %>% select(., -c('pos2', 'ID', 'mate'))
    
    # sorting
    vcf_filter <- vcf 
    vcf_filter$POS1 <- as.integer(vcf_filter$POS1)
    vcf_filter$POS2 <- as.integer(vcf_filter$POS2)
    vcf_filter$chr1_index <- sapply(vcf_filter$CHR1, chr_to_int)
    vcf_filter$chr2_index <- sapply(vcf_filter$CHR2, chr_to_int)
    
    wrong_order1 <- vcf_filter %>% 
      dplyr::filter(CHR1 != CHR2) %>% 
      dplyr::filter(as.numeric(chr1_index) > as.numeric(chr2_index))
    wrong_order2 <- vcf_filter %>% 
      dplyr::filter(CHR1 == CHR2) %>% 
      dplyr::filter(POS1 > POS2)
    wrong <- rbind(wrong_order1, wrong_order2)
    if (nrow(wrong) != 0){
      sorted <- wrong
      sorted$CHR1 <- wrong$CHR2
      sorted$POS1 <- wrong$POS2
      sorted$CHR2 <- wrong$CHR1
      sorted$POS2 <- wrong$POS1
      sorted$dir1 <- wrong$dir2
      sorted$dir2 <- wrong$dir1
      
      correct_order <- vcf_filter %>% 
        dplyr::filter((CHR1 != CHR2 & (as.numeric(chr1_index) <= as.numeric(chr2_index))) | 
                        ((CHR1 == CHR2) & (POS1 <= POS2))) 
      vcf_sorted <- rbind(correct_order, sorted) %>% 
        select(-c('chr1_index', 'chr2_index'))
    } else {
      vcf_sorted <- vcf_filter %>% 
        select(-c('chr1_index', 'chr2_index'))
    }
    
    # remove duplicates 
    vcf_rmdup <- vcf_sorted[!duplicated(vcf_sorted), ]
    vcf_rmdup <- vcf_rmdup[order(vcf_rmdup$CHR1, vcf_rmdup$POS1), ]
    if (nrow(vcf_rmdup) != 1){
      vcf_rmdup$dup <- F
      for (i in 1:(nrow(vcf_rmdup)-1)){
        if (vcf_rmdup$CHR1[i] == vcf_rmdup$CHR1[i+1] & 
            abs(vcf_rmdup$POS1[i] - vcf_rmdup$POS1[i+1]) <= 5 &
            vcf_rmdup$CHR2[i] == vcf_rmdup$CHR2[i+1] & 
            abs(vcf_rmdup$POS2[i] - vcf_rmdup$POS2[i+1]) <= 5 &
            vcf_rmdup$dir1[i] == vcf_rmdup$dir1[i+1] &
            vcf_rmdup$dir2[i] == vcf_rmdup$dir2[i+1]){
          vcf_rmdup$dup[i+1] <- T
        }
      }
      vcf_rmdup <- vcf_rmdup %>% filter(dup == F) %>% select(.,-c('dup'))
    }
    
    # sv type
    vcf_rmdup$SVTYPE <- '.'
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 != vcf_rmdup$CHR2, 'TRA', vcf_rmdup$SVTYPE)
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 == vcf_rmdup$CHR2 & vcf_rmdup$dir1 == vcf_rmdup$dir2, 'INV', vcf_rmdup$SVTYPE)
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 == vcf_rmdup$CHR2 & vcf_rmdup$dir1 == 5 & vcf_rmdup$dir2 == 3, 'DUP', vcf_rmdup$SVTYPE)
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 == vcf_rmdup$CHR2 & vcf_rmdup$dir1 == 3 & vcf_rmdup$dir2 == 5, 'DEL', vcf_rmdup$SVTYPE)
    
    # terminal
    vcf_form <- unite(vcf_rmdup, 'CT', c('dir1', 'dir2'), sep = 'to')
    vcf_form <- select(vcf_form, c('CHR1', 'POS1', 'CHR2', 'POS2', 'SVTYPE', 'CT'))
    
    return(vcf_form)
  }
}
