# GREMLIN: Genomic REarrangements by Machine Learning-based INspection

Automatically refine somatic genomic rearrangements from whole-genome sequences of tumor and matched normal tissues. GREMLIN takes a SV call set as an input and extracts features from the whole-genome bam (or cram). Then, it scores each of the calls and outputs a refined SV list.

GREMLIN was trained and verified using >200k SVs from ~1,800 cancer whole-genomes obtained from the [PCAWG](https://www.nature.com/articles/s41586-019-1913-9) and [Lee et al.](https://www.sciencedirect.com/science/article/pii/S0092867419305112) You can simply apply GREMLIN optimized for the PCAWG dataset or retrain the model with the curated SV calls from a small fraction of samples from your cohort.

## Installation
```
git clone https://github.com/phansol/gremlin/gremlin.git
cd gremlin

Rscript requirements.R
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
* ``-v`` Structural variation call set (tab-delimited; with a header CHR1 POS1 CHR2 POS2 SVTYPE CT<br>(SVTYPE: DEL|DUP|INV|TRA, CT: 3to3|3to5|5to3|5to5)
* ``-n`` Normal bam (or cram)
* ``-t`` Tumor bam (or cram)
* ``-r`` Reference fasta (index should be [given_fasta].fai)
* ``-i`` Sample ID [default: basename of input vcf]
* ``-o`` Output directory [default: directory of input vcf]
* ``-g`` Reference version (19|38) [default: 19]
* ``-c`` Tumor cell fraction [default: 0.5]
* ``-p`` Tumor genome ploidy [default: 2]
* ``-w`` Whole-genome duplication status (wgd|no_wgd) [default: no_wgd]
* ``-y`` Tumor tissue (Biliary|Bladder|Bone_SoftTissue|Breast|Cervix|CNS|Colon_Rectum|Esophagus|Head_Neck|<br>Hematologic|Kideny|Liver|Lung|Ovary|Pancreas|Prostate|Skin|Stomach|Thyroid|Uterus) [default: Biliary] 

## Output
* ``*.feature.dummies.pon.score``: scored SV call set

* ``*.sv.gremlin.tier1.vcf``: tier1 SV calls predicted to be true somatic mutations by GREMLIN

* ``*.sv.gremlin.tier2.vcf``: tier2 SV calls refined with more lenient filtering threshold than tier1

## Best practice
|Step|Description|
|:--:|--|
|1|*(Optional)* [Quality control of input sequences](#quality-control-of-input-sequences)|
|2|[Preprocessing of input SV call sets](#preprocessing-of-input-sv-call-sets)|
|3|[Applying GREMLIN](#usage)|
|4|*(Optional)* [Adjusting classification threshold](#adjusting-classification-threshold)<br>[Retraining GREMLIN](#retraining-gremlin)<br>[Additional filtering using normal panels of your cohort](#additional-filtering-using-normal-panels-of-your-cohort)|

## Quality control of input sequences
Before checking the quality of input sequencing data, install required packages using `Rscript requirements.qc.R`

#### 1. Flag for short inversion artifacts
Short inversion artifacts are a main source of false-positive SV calls, commonly seen in whole-genome sequences of low-quality genomic DNA. Thus, we recommend checking the fraction of short inversions in your sequencing data before applying GREMLIN. 

The following command will estimate the fraction of short inversions among total read pairs using samtools. If your data has an exceptionally high fraction of short inversions, you will get a fail flag, and the refined list (GREMLIN’s output) may include many short inversion errors.
```
Usage: 1_quality_check/samtools_short_inv.sh [TUMOR_BAM/CRAM] [THREADS]

Output: [TUMOR_BAM/CRAM].shinv.pass or [TUMOR_BAM/CRAM].shinv.fail
```

#### 2. Flag for variable sequencing coverage

```
Usage for bam: 1_quality_check/indexcov_read_depth.bam.sh [TUMOR_BAM] [NORMAL_BAM] [OUTPUT_DIRECTORY] [REFERENCE_BUILD] [REFERENCE_FASTA_INDEX]
Usage for cram: 1_quality_check/indexcov_read_depth.cram.sh [TUMOR_CRAM] [NORMAL_CRAM] [OUTPUT_DIRECTORY] [REFERENCE_BUILD] [REFERENCE_FASTA_INDEX]

Output: [TUMOR_BAM/CRAM].depth_ratio.png
        [TUMOR_BAM/CRAM].depth.pass or [TUMOR_BAM/CRAM].depth.fail
```
* ``REFERENCE_BUILD``: reference genome version (19|38)
* ``REFERENCE_FASTA_INDEX``: /path/to/reference.fasta.fai


## Preprocessing of input SV call sets
Formatting SV call sets

If you called SVs using DELLY, SvABA, BRASS, or dRanger, run the following command.
```
Usage: Rscript 2_preprocessing/vcf_formatting.R [VCF] [CALLER] [REFERENCE_FASTA_INDEX] [OUTPUT_DIRECTORY]

Output: [OUTPUT_DIRECTORY]/[VCF].sort
```
* ``VCF``: 
* ``CALLER``: 

Otherwise, transform your SV call set into the ...

tab-separated 

## Adjusting classification threshold
You can adjust filtering threshold (default is 0.89 for tier1 and 0.57 for tier2)
```
Usage: Rscript 5_postprocessing/optional_adjusting_classification_threshold.R [OUTPUT] [THRESHOLD]

Output: [OUTPUT].sv.gremlin.[THRESHOLD].vcf
```
* ``OUTPUT``: feature.dummies.pon.score
* ``THRESHOLD``: classification threshold between 0 and 1 

## Retraining GREMLIN
Before retraining the model, install required packages using `Rscript requirements.rt.R`

#### 1. Re-training with your data
```
Usage: Rscript 5_postprocessing/optional_re_training.R [NEW_DATASET] [OUTPUT_DIRECTORY] [PREFIX] [PERCENT]

Output: [OUTPUT_DIRECTORY]/[PREFIX]_gbm.fit.rds
```
* ``NEW_DATASET``: same format with ``feature.dummies.pon.score`` with an additional column (true_label = T/F)
* ``PERCENT``: (optional) the percent of the training set will be used for re-training; value between 0 and 100
	                  160 training samples will be used in default
	       
#### 2. Applying the re-trained model to your data
```
Usage: Rscript 5_postprocessing/optional_apply_re_trained.R [OUTPUT] [re_trained_gbm.fit.rds] [THRESHOLD]

Output: [OUTPUT].re_trained
	[OUTPUT].re_trained.[THRESHOLD].vcf (if threshold is given)
```
* ``OUTPUT``: feature.dummies.pon.score
* ``THRESHOLD``: (optional) used for filtering vcf; value between 0 and 1

## Additional filtering using normal panels of your cohort
Cohort-specific panel of normals 

## License
