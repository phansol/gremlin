#' args1 = bed1.seq
#' args2 = bed2.seq
#' args3 = vcf

args <- commandArgs(trailingOnly = T)
options(warn = -1)

suppressMessages(library(readr))
suppressMessages(library(stringr))

log_file <- paste0(dirname(args[3]), '/log/', basename(args[3]) %>% gsub('.vcf.sort.SVinfo.fi', '', .), '.s3.success')
out_file <- paste0(args[3], '.gc_mask')

if (!(file.exists(log_file) & file.exists(out_file))){
  
  bed1 <- read_tsv(args[1], col_names = F)
  bed2 <- read_tsv(args[2], col_names = F)
  vcf <- read_tsv(args[3], col_types = cols(.default = 'c'))
  
  if (nrow(bed1) == nrow(bed2) & nrow(bed2) == nrow(vcf)){
    for (j in 1:nrow(bed1)){
      # BP1
      N1 <- str_count(bed1[j, 2], 'N') + str_count(bed1[j, 2], 'n')
      A1 <- str_count(bed1[j, 2], 'A')
      T1 <- str_count(bed1[j, 2], 'T')
      G1 <- str_count(bed1[j, 2], 'G')
      C1 <- str_count(bed1[j, 2], 'C')
      a1 <- str_count(bed1[j, 2], 'a')
      t1 <- str_count(bed1[j, 2], 't')
      g1 <- str_count(bed1[j, 2], 'g')
      c1 <- str_count(bed1[j, 2], 'c')
      
      if (nchar(bed1[j, 2]) != N1){
        gc1 <- (G1 + C1 + g1 + c1)/(nchar(bed1[j, 2]) - N1)
      } else {
        gc1 <- 0
      }
      
      sm1 <- (N1 + a1 + t1 + g1 + c1)/nchar(bed1[j, 2]) 
      vcf$gc_content1[j] <- round(gc1, 4)
      vcf$soft_masked1[j] <- round(sm1, 4)
      
      # BP2
      N2 <- str_count(bed2[j, 2], 'N') + str_count(bed2[j, 2], 'n')
      A2 <- str_count(bed2[j, 2], 'A')
      T2 <- str_count(bed2[j, 2], 'T')
      G2 <- str_count(bed2[j, 2], 'G')
      C2 <- str_count(bed2[j, 2], 'C')
      a2 <- str_count(bed2[j, 2], 'a')
      t2 <- str_count(bed2[j, 2], 't')
      g2 <- str_count(bed2[j, 2], 'g')
      c2 <- str_count(bed2[j, 2], 'c')
      
      if (nchar(bed2[j, 2]) != N2){
        gc2 <- (G2 + C2 + g2 + c2)/(nchar(bed2[j, 2]) - N2)
      } else {
        gc2 <- 0
      }
      sm2 <- (N2 + a2 + t2 + g2 + c2)/nchar(bed2[j, 2]) 
      vcf$gc_content2[j] <- round(gc2, 4)
      vcf$soft_masked2[j] <- round(sm2, 4)
    }
  } else {
    print('ERROR: the number of lines in bed and vcf are not the same!')
  }
  write_tsv(vcf, out_file)
  
}
