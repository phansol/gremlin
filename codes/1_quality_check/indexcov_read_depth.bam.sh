#!/bin/bash
tbam=$1
nbam=$2
outdir=$3
refver=$4 # either 19 or 38
reffai=$5

goleft=/home/users/hspark/Tools/goleft

set -e

# run indexcov
$goleft indexcov -d $outdir/$(basename $tbam).tmp $tbam
$goleft indexcov -d $outdir/$(basename $nbam).tmp $nbam

# post processing
Rscript "$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/indexcov_read_depth_postprocessing.R" \
	$outdir/$(basename $tbam).tmp/$(basename $tbam).tmp-indexcov.bed.gz \
	$outdir/$(basename $nbam).tmp/$(basename $nbam).tmp-indexcov.bed.gz \
	$refver \
	$reffai

mv $outdir/$(basename $tbam).tmp/$(basename $tbam).tmp-indexcov.bed.gz $outdir/$(basename $tbam).indexcov.bed.gz
mv $outdir/$(basename $nbam).tmp/$(basename $nbam).tmp-indexcov.bed.gz $outdir/$(basename $nbam).indexcov.bed.gz
rm -r $outdir/$(basename $tbam).tmp $outdir/$(basename $nbam).tmp
