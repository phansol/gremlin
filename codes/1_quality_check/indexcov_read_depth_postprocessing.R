#' variable coverage flag
#' args1 = tumour.indexcov.bed.gz
#' args2 = normal.indexcov.bed.gz
#' args3 = reference version (19 or 38)
#' args4 = reference.fasta.fai

suppressMessages(library(tidyverse))
suppressMessages(library(GWASTools))
suppressMessages(library(rpart))

options(warn = -1)

args <- commandArgs(trailingOnly = T)
  
message('# estimating the variance of read depths')

# inputs
tcov <- read_tsv(args[1], col_types = cols(`#chrom` = 'c'))
ncov <- read_tsv(args[2], col_types = cols(`#chrom` = 'c'))

if (args[3] == '19'){
  data(centromeres.hg19)
  centromere <- centromeres.hg19
} else if (args[3] == '38'){
  data(centromeres.hg38)
  centromere <- centromeres.hg38
}

fai <- read_tsv(args[4], col_names = c('chr', 'size', 'x', 'y', 'z'), col_types = cols())

cov_trim <- function(cov){
  
  # match 'chr' prefix
  if (grepl('chr', cov$`#chrom`[1])){
    centromere$chrom <- centromere$chrom %>% gsub('chr', '', .) %>% paste0('chr', .)
    fai$chr <- fai$chr %>% gsub('chr', '', .) %>% paste0('chr', .)
  } else {
    centromere$chrom <- centromere$chrom %>% gsub('chr', '', .) 
    fai$chr <- fai$chr %>% gsub('chr', '', .)
  }
  
  # primary chromosomes
  cov <- cov %>% filter(`#chrom` %in% c(c(1:22, 'X'), paste0('chr', c(1:22, 'X')))) 
  
  # exclude centomeric & telomeric regions
  for (i in 1:23){
    cov <- cov %>%
      filter(!(`#chrom` == centromere$chrom[i] & 
                 ((start > centromere$left.base[i] & start < centromere$right.base[i]) | (end > centromere$left.base[i] & end < centromere$right.base[i]))))
    cov <- cov %>%
      filter(!(`#chrom` == fai$chr[i] & ((end > fai$size[i] - 100000) | (start < 100000))))
    
    # acrocentric or telocentric chromosomes (chr13, 14, 15, 21, 22)
    if (i %in% c(13, 14, 15, 21, 22)){
      cov <- cov %>%
        filter(!(`#chrom` == centromere$chrom[i] & end < centromere$right.base[i]))
    }
  }
  cov
}

tcov <- cov_trim(tcov)
ncov <- cov_trim(ncov)


# depth ratio
colnames(tcov)[4] <- 'tumor'
colnames(ncov)[4] <- 'normal'
cov <- merge(tcov, ncov, by = c('#chrom', 'start', 'end')) %>% as.data.frame()
cov$cn_estimate <- cov$tumor/cov$normal*2
cov$`#chrom` <- cov$`#chrom` %>% gsub('chr', '', .) %>% gsub('X', 23, .) %>% as.numeric()
cov <- cov[order(cov$`#chrom`, cov$start), ]
cov <- cov[complete.cases(cov), ]


# remove extreme cn cases
m <- quantile(cov$cn_estimate, 0.005)
M <- quantile(cov$cn_estimate, 0.995)
cov <- cov %>% filter(cn_estimate >= m & cn_estimate <= M)


# cnv detection using recursive partitioning
cov$index <- 1:nrow(cov)
cov$cn_smoothened <- 0

sample <- basename(args[1]) %>% gsub('.tmp-indexcov.bed.gz', '', .)
sample <- ifelse(nchar(sample) <= 25, sample,
                 substring(sample, seq(1, nchar(sample)-1, 25), unique(c(seq(25, nchar(sample), 25), nchar(sample)))) %>% paste0(., collapse = '\n'))

png(dirname(args[1]) %>% gsub('tmp$', 'depth_ratio.png', .), height = 2970, width = 2100)
par(mfrow = c(8, 3), mar = c(8, 8, 4, 4), mgp = c(6, 2, 0))

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, sample, cex = 5)

for (i in 1:23){
  cov_i <- cov %>% filter(`#chrom` == i)
  tree <- rpart(cn_estimate ~ start, data = cov_i, minsplit = 1, minbucket = 60) # minbucket = 60 corressponds to 1Mb
  cov_i$cn_smoothened <- predict(tree, data.frame(start = cov_i$start))
  cov$cn_smoothened[cov_i$index] <- cov_i$cn_smoothened
  plot(cov_i$start, cov_i$cn_estimate, ylim = c(max(min(cov$cn_estimate) -0.2, 0), max(cov$cn_estimate) + 0.2), cex.axis = 2, cex.lab = 3, 
       xlab = paste0('chromosome ', i %>% gsub(23, 'X', .)), ylab = 'Read depth ratio x 2')
  lines(cov_i$start, cov_i$cn_smoothened, col = 'red', lwd = 5)
}

p <- dev.off()
message(paste0('# see ', dirname(args[1]) %>% gsub('tmp$', 'depth_ratio.png', .)))

# variance of read depths (kb bin)
v <- sum((cov$cn_smoothened - cov$cn_estimate)^2)/nrow(cov)

if (v > 0.05){
  writeLines(paste0('# variance of read depths: ', v), dirname(args[1]) %>% gsub('tmp$', 'depth.fail', .))
  message(paste0('# See ', dirname(args[1]) %>% gsub('tmp$', 'depth.fail', .), ' for the estimate of variance of read depths'))
  message(paste0('# WARNING: ', dirname(args[1]) %>% gsub('.tmp$', '', .), ' is suspected of highly variable read depths'))
} else {
  message(paste0('# See ', dirname(args[1]) %>% gsub('tmp$', 'depth.pass', .), ' for the estimate of variance of read depths'))
  writeLines(paste0('# variance of read depths: ', v), dirname(args[1]) %>% gsub('tmp$', 'depth.pass', .))
}

message('# Done')
