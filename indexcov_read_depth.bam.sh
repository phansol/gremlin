#!/bin/bash
tbam=$1
nbam=$2
outdir=$3
refver=$4 # either 19 or 38
reffai=$5

# reffai
# hg19: /home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai
# hg38: /home/users/pjh/References/reference_genome/GRCh38/GCA_for_alignment_pipelines/no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai


goleft=/home/users/hspark/Tools/goleft

set -e

# run indexcov
$goleft indexcov -d $outdir/$(basename $tbam).tmp $tbam
$goleft indexcov -d $outdir/$(basename $nbam).tmp $nbam

# post processing
Rscript /home/users/hspark/Scripts/goleft/indexcov_read_depth_postprocessing.R \
	$outdir/$(basename $tbam).tmp/$(basename $tbam).tmp-indexcov.bed.gz \
	$outdir/$(basename $nbam).tmp/$(basename $nbam).tmp-indexcov.bed.gz \
	$refver \
	$reffai

mv $outdir/$(basename $tbam).tmp/$(basename $tbam).tmp-indexcov.bed.gz $outdir/$(basename $tbam).indexcov.bed.gz
mv $outdir/$(basename $nbam).tmp/$(basename $nbam).tmp-indexcov.bed.gz $outdir/$(basename $nbam).indexcov.bed.gz
rm -r $outdir/$(basename $tbam).tmp $outdir/$(basename $nbam).tmp
