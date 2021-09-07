#' args1 : vcf
#' args2 : reference fasta.fai

args <- commandArgs(trailingOnly = T)
options(warn = -1)

suppressMessages(library(tidyverse))

log_file <- paste0(dirname(args[1]), '/log/', basename(args[2]) %>% gsub('.vcf.sort.SVinfo.fi', '', .), '.s3.success')

if (!(file.exists(log_file) & file.exists(paste0(args[1], '.gc_mask')))){

  fai <- read_tsv(args[2], col_names = F) # chr, size
  vcf <- read_tsv(args[1], col_types = cols(.default = 'c'))
  vcf$index <- 1:nrow(vcf)
  
  bed1 <- select(vcf, c('index', 're_chr1', 're_pos1')) 
  bed1$start <- as.integer(vcf$re_pos1) - 50
  bed1$end <- as.integer(vcf$re_pos1) + 51
  colnames(bed1)[2] <- 'chr'
  
  bed2 <- select(vcf, c('index', 're_chr2', 're_pos2')) 
  bed2$start <- as.integer(vcf$re_pos2) - 50
  bed2$end <- as.integer(vcf$re_pos2) + 51
  colnames(bed2)[2] <- 'chr'
  
  correct_position <- function(bed, fai){
    for (i in 1:nrow(bed)){
      bed$start <- as.integer(bed$start)
      bed$end <- as.integer(bed$end)
      if (bed$start[i] < 0){
        bed$end[i] <- bed$end[i] - bed$start[i] + 1
        bed$start[i] <- 1
      } else if (bed$end[i] > fai$X2[which(bed$chr[i] == fai$X1)]){
        bed$start[i] <- bed$start[i] - (bed$end[i] - fai$X2[which(bed$chr[i] == fai$X1)])
        bed$end[i] <- fai$X2[which(bed$chr[i] == fai$X1)] 
      }
    }
    bed <- bed[order(bed$index), ]
    bed <- select(bed, c('chr', 'start', 'end'))
    return(bed)
  }
  
  bed1 <- correct_position(bed1, fai)
  bed2 <- correct_position(bed2, fai)
  
  bed1$start <- as.integer(bed1$start)
  bed1$end <- as.integer(bed1$end)
  bed2$start <- as.integer(bed2$start)
  bed2$end <- as.integer(bed2$end)
  
  write_tsv(bed1, paste0(args[1], '.bed1'), col_names = F)
  write_tsv(bed2, paste0(args[1], '.bed2'), col_names = F)
  
}

# with output files, do the next command
# bedtools getfasta -fi [fasta] -bed [bed_file] -tab -fo [output_file_name]
