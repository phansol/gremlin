#' formatting vcf
#' args1 = vcf
#' args2 = sv caller [delly/svaba/brass/dranger/lumpy/manta]
#' args3 = reference.fa.fai
#' args4 = output directory [default = vcf directory]

suppressMessages(library(tidyverse))
options(warn = -1)

args <- commandArgs(trailingOnly = F)
srcdir <- args[which(grepl('--file=', args))] %>% gsub('--file=', '', .) %>% dirname()
args <- args[(which(grepl('--args', args))+1):length(args)]

if (args[2] == 'delly'){
  source(paste0(srcdir, '/delly.R'))
} else if (args[2] == 'svaba'){
  source(paste0(srcdir, '/svaba.R'))
} else if (args[2] == 'brass'){
  source(paste0(srcdir, '/brass.R'))
} else if (args[2] == 'dranger'){
  source(paste0(srcdir, '/dranger.R'))
} else if (args[2] == 'lumpy'){
  source(paste0(srcdir, '/lumpy.R'))
} else if (args[2] == 'manta'){
  source(paste0(srcdir, '/manta.R'))
} else {
  message('# WARNING: wrong input')
  quit()
}

vcf_form <- vcf_formatting(vcf_path = args[1], ref.fai_path = args[3])

outdir <- ifelse(length(args) == 4, args[4], dirname(args[1]))
system(paste0('mkdir -p ', outdir))

if (vcf_form == 'empty'){
  message('# Input vcf is empty')
  message('# Job finished')
} else {
  outfile <- paste0(outdir, '/', basename(args[1]) %>% gsub('.gz', '', .), '.sort')
  write_tsv(vcf_form, outfile)
  message(paste0('# See ', outfile))
  message('# Job finished')
}
