#!/bin/bash
tcrai=$1
ncrai=$2
outdir=$3
refver=$4 # either 19 or 38
reffai=$5

goleft=/home/users/hspark/Tools/goleft

set -e

# run indexcov
$goleft indexcov --extranormalize -d $outdir/$(basename $tcrai).tmp --fai $reffai $tcrai
$goleft indexcov --extranormalize -d $outdir/$(basename $ncrai).tmp --fai $reffai $ncrai

# post processing
Rscript /home/users/hspark/Scripts/goleft/indexcov_read_depth_postprocessing.R \
	$outdir/$(basename $tcrai).tmp/$(basename $tcrai).tmp-indexcov.bed.gz \
	$outdir/$(basename $ncrai).tmp/$(basename $ncrai).tmp-indexcov.bed.gz \
	$refver \
	$reffai

mv $outdir/$(basename $tcrai).tmp/$(basename $tcrai).tmp-indexcov.bed.gz $outdir/$(basename $tcrai).indexcov.bed.gz
mv $outdir/$(basename $ncrai).tmp/$(basename $ncrai).tmp-indexcov.bed.gz $outdir/$(basename $ncrai).indexcov.bed.gz
rm -r $outdir/$(basename $tcrai).tmp $outdir/$(basename $ncrai).tmp
