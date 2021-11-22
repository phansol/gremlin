# LUMPY 

chr_to_int <- function(x){
  y <- x %>% gsub('chr', '', .) %>% gsub('X', 23, .) %>% gsub('Y', 24, .) %>% as.numeric()
  return(y)
}

chr_size <- function(chr, fai){
  s <- fai$X2[which(fai$X1 == chr)]
  return(s)
}

vcf_formatting <- function(vcf_path, ref.fai_path){
  primary_chr <- c(1:22, 'X', 'Y', paste0('chr', c(1:22, 'X', 'Y')))
  
  # input
  vcf <- read_tsv(vcf_path, col_names = c('CHR1', 'POS1', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR'), comment = '#',
                  col_types = cols('CHR1' = 'c')) %>%
    separate(INFO, into = c('a', 'b'), sep = 'STRANDS=', remove = F) %>%
    separate(b, into = c('CT', 'c'), sep = ':')
    
  vcf$dir1 <- ifelse(grepl('^\\+', vcf$CT), 3, 5)
  vcf$dir2 <- ifelse(grepl('\\+$', vcf$CT), 3, 5)

  vcf_single <- vcf %>% 
    filter(ALT %in% c('<DEL>', '<DUP>', '<INV>')) %>% 
    separate(INFO, into = c('d', 'e'), sep = ';END=') %>%
    separate(e, into = c('POS2', 'f'), sep = ';') %>%
    mutate(CHR2 = CHR1) %>%
    select(., c('CHR1', 'POS1', 'CHR2', 'POS2', 'dir1', 'dir2')) %>%
    filter(CHR1 %in% primary_chr & CHR2 %in% primary_chr)
  
  vcf_pair <- vcf %>% 
    filter(! ALT %in% c('<DEL>', '<DUP>', '<INV>')) %>% 
    separate(ALT, into = c('chr2', 'pos2'), sep = ':') %>%
    separate(chr2, into = c('a', 'CHR2'), sep = '\\[|\\]') %>%
    separate(pos2, into = c('POS2', 'b'), sep = '\\[|\\]') %>%
    select(c('CHR1', 'POS1', 'CHR2', 'POS2', 'dir1', 'dir2')) %>%
    filter(CHR1 %in% primary_chr & CHR2 %in% primary_chr)
  
  vcf <- rbind(vcf_single, vcf_pair)
  
  if (nrow(vcf) == 0){
    return('empty')
  } else {
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
    
    # sv type
    vcf_rmdup$SVTYPE <- '.'
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 != vcf_rmdup$CHR2, 'TRA', vcf_rmdup$SVTYPE)
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 == vcf_rmdup$CHR2 & vcf_rmdup$dir1 == vcf_rmdup$dir2, 'INV', vcf_rmdup$SVTYPE)
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 == vcf_rmdup$CHR2 & vcf_rmdup$dir1 == 5 & vcf_rmdup$dir2 == 3, 'DUP', vcf_rmdup$SVTYPE)
    vcf_rmdup$SVTYPE <- ifelse(vcf_rmdup$CHR1 == vcf_rmdup$CHR2 & vcf_rmdup$dir1 == 3 & vcf_rmdup$dir2 == 5, 'DEL', vcf_rmdup$SVTYPE)
    
    # terminal
    vcf_rmdup <- vcf_rmdup %>%
      unite(., 'CT', c('dir1', 'dir2'), sep = 'to') %>%
      select(., c('CHR1', 'POS1', 'CHR2', 'POS2', 'SVTYPE', 'CT'))
    
    # filter if a breakpoint locates at extreme end of the chromosome
    fai <- read_tsv(ref.fai_path, col_names = F)
    
    vcf_rmdup$chr1_size <- sapply(vcf_rmdup$CHR1, function(x) chr_size(x, fai))
    vcf_rmdup$chr2_size <- sapply(vcf_rmdup$CHR2, function(x) chr_size(x, fai))
    
    vcf_form <- vcf_rmdup %>%
      filter(as.numeric(POS1) > 1000 & as.numeric(POS2) > 1000) %>%
      filter(as.numeric(POS1) < chr1_size - 1000 & as.numeric(POS2) < chr2_size - 1000) %>%
      select(., -c('chr1_size', 'chr2_size'))
    
    return(vcf_form)
  }
} 
