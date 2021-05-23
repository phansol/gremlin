# GREMLIN: Genomic REarrangements by Machine Learning-based INspection

Automatically refine somatic genomic rearrangements from whole-genome sequences of tumor and matched normal tissues. GREMLIN takes a SV call set as an input and extracts features from the whole-genome bam (or cram). Then, it scores each of the calls and outputs a refined SV list.

GREMLIN was trained and verified using >200k SVs from ~1,800 cancer whole-genomes obtained from the PCAWG and Lee et al. You can simply apply GREMLIN optimized for the PCAWG dataset or retrain the model with the curated SV calls from a small fraction of samples from your cohort.

Table of contents
=================

  * [Installation]
  * [Usage]
  * [Output file description]
  * [Adjusting threshold]
  * [Retraining GREMLIN]

## Installation

git clone https://github.com/phansol/gremlin/gremlin.git
cd gremlin

Rscript install_R_requirements.R
pip install -r requirements.txt

## Usage
```
gremlin [-h] [-v VCF] [-n NORMAL_BAM] [-t TUMOR_BAM] [-r REFERENCE_FASTA]
        [-i SAMPLE_ID] [-o OUTPUT_DIRECTORY] [-g REFERENCE_VERSION]
        [-c TUMOR_CELL_FRACTION] [-p TUMOR_PLOIDY] [-w WGD_STATUS] [-y TUMOR_TISSUE]
```
Required arguments:
* ``-v`` List of structural variations in the form of CHR1/POS1/CHR2/POS2/SVTYPE/CT
         (SVTYPE: DEL|DUP|INV|TRA, CT: 3to3|3to5|5to3|5to5)
* ``-n`` Normal bam (or cram)
* ``-t`` Tumor bam (or cram)
* ``-r`` Reference fasta (index should be [given_fasta].fai)
* ``-i`` Sample ID [default: basename of input vcf]
* ``-o`` Output directory [default: directory of input vcf]
* ``-g`` Reference version (19|38) [default: 19]
* ``-c`` Tumor cell fraction [default: 0.5]
* ``-p`` Tumor genome ploidy [default: 2]
* ``-w`` Whole-genome duplication status (wgd|no_wgd) [default: no_wgd]
* ``-y`` Tumor tissue (Biliary|Bladder|Bone_SoftTissue|Breast|Cervix|CNS|Colon_Rectum|Esophagus|
         Head_Neck|Hematologic|Kideny|Liver|Lung|Ovary|Pancreas|Prostate|Skin|Stomach|Thyroid|Uterus) 
         [default: Biliary] 

## Output
##### ``*.feature.dummies.pon.score``
feature annotated vcf

##### ``*.sv.gremlin.tier1.vcf``
GREMLIN 

##### ``*.sv.gremlin.tier2.vcf``


## Test run
```
gremlin -v test/sv.vcf.sort -n test/normal.chr21.bam -t test/tumor.chr21.bam -r test/grch37.chr21.fasta
```

## Adjusting threshold
change threshold

## Retraining GREMLIN

## License
