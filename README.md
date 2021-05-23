# GREMLIN: Genomic REarrangements by Machine Learning-based INspection

Automatically refine somatic genomic rearrangements from whole-genome sequences of tumor and matched normal tissues. GREMLIN takes a SV call set as an input and extracts features from the whole-genome bam (or cram). Then, it scores each of the calls and outputs a refined SV list.

GREMLIN was trained and verified using >200k SVs from ~1,800 cancer whole-genomes obtained from the PCAWG and Lee et al. You can simply apply GREMLIN optimized for the PCAWG dataset or retrain the model with the curated SV calls from a small fraction of samples from your cohort.

## Table of contents
=================

  * [Installation]
  * [Usage]
  * [Output file description]
  * [Adjusting threshold]
  * [Retraining GREMLIN]

## Installation
```
git clone https://github.com/phansol/gremlin/gremlin.git
cd gremlin

Rscript install_R_requirements.R
pip install -r requirements.txt
```

## Test run
```
gremlin -v test/sv.vcf.sort -n test/normal.chr21.bam -t test/tumor.chr21.bam -r test/grch37.chr21.fasta
```

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
* ``-y`` Tumor tissue (Biliary|Bladder|Bone_SoftTissue|Breast|Cervix|CNS|Colon_Rectum|Esophagus|Head_Neck|Hematologic|Kideny|Liver|Lung|Ovary|Pancreas|Prostate|Skin|Stomach|Thyroid|Uterus) [default: Biliary] 

## Output
##### ``*.feature.dummies.pon.score``
scored SV call set

##### ``*.sv.gremlin.tier1.vcf``
tier1 SV calls predicted to be true somatic mutations by GREMLIN

##### ``*.sv.gremlin.tier2.vcf``
tier2 SV calls refined with more lenient filtering threshold than tier1

## Quality control of input sequencing data
### 1. Flag for short inversion artifacts
Short inversion artifacts are a main source of false-positive SV calls, commonly seen in whole-genome sequences of low-quality genomic DNA. Thus, we recommend checking the fraction of short inversions in your sequencing data before applying GREMLIN. 

The following command will estimate the fraction of short inversions among total read pairs using samtools. If your data has an exceptionally high fraction of short inversions, you will get a fail flag, and the refined list (GREMLIN’s output) may include many short inversion errors.
```
Usage: 1_quality_check/samtools_short_inv.sh [TUMOR_BAM] [THREADS]

Output: [TUMOR_BAM].shinv.pass or [TUMOR_BAM].shinv.fail
```

### 2. Flag for variable sequencing coverage

```
Usage for bam: 1_quality_check/indexcov_read_depth.bam.sh [TUMOR_BAM] [NORMAL_BAM] [OUTPUT_DIRECTORY] [REFERENCE_BUILD] [REFERENCE_FASTA_INDEX]
Usage for cram: 1_quality_check/indexcov_read_depth.cram.sh [TUMOR_CRAM] [NORMAL_CRAM] [OUTPUT_DIRECTORY] [REFERENCE_BUILD] [REFERENCE_FASTA_INDEX]

Output: [TUMOR_BAM/CRAM].depth_ratio.png
        [TUMOR_BAM/CRAM]..depth.pass or tumor.bam.depth.fail
```

## Formatting SV call sets

## Cohort-specific panel of normals 


## Adjusting classification threshold
You can adjust filtering threshold (default is 0.89 for tier1 and 0.57 for tier2)
```
Usage: Rscript 5_postprocessing/optional_adjusting_classification_threshold.R [feature.dummies.pon.score] [threshold]

[feature.dummies.pon.score]: 
[threshold]: 

Output: [feature.dummies.pon.score].sv.gremlin.[threshold].vcf
```


## Retraining GREMLIN
### 1. Re-training with your data
```
Usage: Rscript 5_postprocessing/optional_re_training.R <new_dataset> [OUTPUT_DIRECTORY] <prefix> <percent (optional)>

new_dataset: same format with *.feature.dummies.pon.score with an additional column (true_label = T/F)
percent (optional): If given a value X (0-100), X% percent of the training set will be used for re-training. 
		   Otherwise, 160 training samples will be used.

Output: [OUTPUT_DIRECTORY]/prefix_gbm.fit.rds
```
### 2. Applying the re-trained model to your data
```
Usage: Rscript 5_postprocessing/optional_apply_re_trained.R [feature.dummies.pon.score] [re_trained_gbm.fit.rds] <threshold (optional)>

threshold (optional): used for filtering vcf, value between (0, 1).

Output: [feature.dummies.pon.score].re_trained
	  *.sv.re_trained.threshold.vcf (if threshold is given)
```

## License
