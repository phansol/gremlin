# functions 

chr_to_int <- function(x){
  y <- x %>% gsub('chr', '', .) %>% gsub('X', 23, .) %>% gsub('Y', 24, .) %>% as.numeric()
  return(y)
}

vcf_formatting <- function(vcf_path, ref.fai_path){
  primary_chr <- c(1:22, 'X', 'Y', paste0('chr', c(1:22, 'X', 'Y')))
  
  # input
  vcf <- read_tsv(vcf_path, col_names = c('CHR1', 'POS1', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL'), comment = '#',
                  col_types = cols('CHR1' = 'c')) %>% 
    separate(INFO, into = c('a', 'b'), sep = 'CHR2=', remove = F) %>%
    separate(b, into = c('CHR2', 'c'), sep = ';') %>%
    separate(INFO, into = c('d', 'e'), sep = ';END=', remove = F) %>%
    separate(e, into = c('POS2', 'f'), sep = ';') %>%
    separate(INFO, into = c('g', 'h'), sep = 'SVTYPE=', remove = F) %>%
    separate(h, into = c('SVTYPE', 'i'), sep = ';') %>%
    separate(INFO, into = c('j', 'k'), sep = 'CT=', remove = F) %>%
    separate(k, into = c('connection_type', 'l'), sep = ';') %>%
    separate(connection_type, into = c('dir1','dir2'), sep = 'to') 
  
  fm <- strsplit(vcf$FORMAT[1], ':')[[1]]
  vcf <- vcf %>%
    separate(NORMAL, into = fm, sep = ':') %>%
    select(c('CHR1', 'POS1', 'CHR2', 'POS2', 'FILTER', 'SVTYPE', 'dir1', 'dir2', 'DV', 'RV')) %>%
    filter(CHR1 %in% primary_chr & CHR2 %in% primary_chr)
  
  vcf$SVTYPE <- ifelse(vcf$SVTYPE == 'BND', 'TRA', vcf$SVTYPE)
  
  # prefilter -- remove when RV+DV >= 1 in normal sample
  vcf_pass <- vcf %>% filter(FILTER == 'PASS')
  vcf_unpass <- vcf %>% filter(FILTER != 'PASS') %>% 
    filter(as.numeric(DV) == 0 & as.numeric(RV) == 0)
  
  # sorting
  vcf_filter <- rbind(vcf_pass, vcf_unpass) %>% 
    select(-c('FILTER', 'DV', 'RV'))
  if (nrow(vcf_filter) == 0){
    return('empty')
  } else {
    vcf_filter$POS1 <- as.integer(vcf_filter$POS1)
    vcf_filter$POS2 <- as.integer(vcf_filter$POS2)
    vcf_filter$chr1_index <- sapply(vcf_filter$CHR1, chr_to_int)
    vcf_filter$chr2_index <- sapply(vcf_filter$CHR2, chr_to_int)
    
    wrong_order1 <- vcf_filter %>% 
      dplyr::filter(SVTYPE == 'TRA') %>% 
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
        dplyr::filter((SVTYPE == 'TRA' & (as.numeric(chr1_index) <= as.numeric(chr2_index))) | 
                        ((CHR1 == CHR2) & (POS1 <= POS2))) 
      vcf_sorted <- rbind(correct_order, sorted) %>% 
        select(-c('chr1_index', 'chr2_index'))
    } else {
      vcf_sorted <- vcf_filter %>% 
        select(-c('chr1_index', 'chr2_index'))
    }
    vcf_rmdup <- vcf_sorted[!duplicated(vcf_sorted), ]
    vcf_form <- unite(vcf_rmdup, 'CT', c('dir1', 'dir2'), sep = 'to')
    
    return(vcf_form)
  }
}

