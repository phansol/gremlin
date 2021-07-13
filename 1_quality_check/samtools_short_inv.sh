#!/bin/bash
bam=$1
threads=$2

# samtools read count
  # -f --require-flags: output reads that fulfill the checked 'bitcode' criteria
  # -F --exclude-flags: exclude reads that match one or more checked 'bitcode' criteria

# total read count
tc=$(samtools view -@ $threads -f 1 -F 3852 -c $bam)
echo $tc
# inversion read count
ttinv=$(samtools view -@ $threads -f 1 -F 3900 $bam | cut -f9 | awk '($0>0 && $0<10000)||($0<0 && $0>-10000)' | wc -l)
hhinv=$(samtools view -@ $threads -f 49 -F 3852 $bam | cut -f9 | awk '($0>0 && $0<10000)||($0<0 && $0>-10000)' | wc -l)
echo $ttinv
echo $hhinv
# inversion fraction
let inv=$ttinv+$hhinv
echo $inv
inv_frac=$(echo $inv $tc | awk '{print $1/$2}')
echo $inv_frac
echo "# short inversion fraction: $inv_frac" > $bam.shinv

if (( $(bc <<< "inv_frac > 0.0015") )); then
    mv $bam.shinv $bam.shinv.fail
    echo "# WARNING: $bam is suspected of many short inversion artifacts"
    echo "# see $bam.shinv.fail for the fraction of short inversions (<10kb)"
else
    mv $bam.shinv $bam.shinv.pass
    echo "# done"
    echo "# see $bam.shinv.pass for the fraction of short inversions (<10kb)"
fi

# SAM FLAG
  # see http://broadinstitute.github.io/picard/explain-flags.html
  # see https://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/ 

  # 1    read paired
  # 2    read mapped in proper pair
  # 4    read unmapped
  # 8    mate unmapped
  # 16   read reverse strand
  # 32   mate reverse strand
  # 64   first in pair
  # 128  second in pair
  # 256  not primary alignment
  # 512  read fails platform/vendor quality checks
  # 1024 read is PCR or optical duplicate
  # 2048 supplementary alignment
