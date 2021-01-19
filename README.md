# GREMLIN: Genomic REarrangements by Machine Learning-based INspection

GREMLIN is a machine learning classifier that can refine somatic structural variation calls. 

Table of contents
=================

  * [Installation]
  * [Usage]
  * [Output file description]

## Installation


## Usage
```
sh gremlin.sh -v <vcf> -o <outdir> -i <sample_id> -n <path/to/normal_bam> -t <path/to/tumor_bam> -r <fasta> -g <19/38> -c <purity> -p <ploidy> -w <wgd_status>
```

Required options:
* ``-v`` vcf, structural variant candidates. [example] 
tab-delimited format
CHR1    POS1    CHR2    POS2    SVTYPE  CT 
* ``-i`` Sample ID
* ``-n`` Normal bam
* ``-t`` Tumor bam
* ``-r`` Reference fasta (the file name of fasta.fai should be [given_fasta_file_name].fai)
* ``-g`` Reference genome version (either 19 or 38)
* ``-c`` cellularity (tumor purity)
* ``-p`` ploidy
* ``-w`` whole-genome duplication status (either wgd or no_wgd)
* ``-y`` tumor histology (Biliary/Bladder/Bone_SoftTissue/Breast/Cervix/CNS/Colon_Rectum/Esophagus/Head_Neck/Hematologic/Kideny/Liver/Lung/Ovary/Pancreas/Prostate/Skin/Stomach/Thyroid/Uterus)
* ``-o`` Output directory [current directory]

## Output
##### ``*.features``
feature annotated vcf

##### ``*.tier1.vcf``
GREMLIN 

##### ``*.tier2.vcf``

## Re-filtering 
change threshold


## License
