[![DOI](https://zenodo.org/badge/330940350.svg)](https://zenodo.org/badge/latestdoi/330940350)
# GREMLIN: Genomic REarrangements by Machine Learning-based INspection

### Machine learning-based detection of structural variations in human cancer genomes

#### Hansol Park, Seongyeol Park, Ji-Hyung Park, Jeonghwan Youk, Kwihoon Kim, Su Yeon Kim, and Young Seok Ju


GREMLIN automatically refines somatic structural variations (SVs) from whole-genome sequences of tumor and matched normal tissues. GREMLIN takes a SV call set as an input and extracts features from the whole-genome bam (or cram). Then, it scores each of the calls and outputs a refined SV list.

GREMLIN was trained and verified using >200k SVs from 1,802 cancer whole-genomes obtained from the [PCAWG](https://www.nature.com/articles/s41586-019-1913-9) and [Lee et al.](https://www.sciencedirect.com/science/article/pii/S0092867419305112) You can simply apply GREMLIN optimized for the PCAWG dataset or retrain the model with the curated SV calls from a small fraction of samples from your cohort.

## System requirements
### Hardware requirements
GREMLIN requires only a standard computer with enough RAM to support the operations. 

### Operating systems
GREMLIN has been tested on *Linux (Ubuntu 18.04 LTS).*

### Software dependencies
To implement GREMLIN, python (v3.6.6), R (v3.6.0), and bedtools (v2.25.0) are required.

GREMLIN depends on the following python and R packages. 
- Python dependencies: *pyranges, pysam*
- R dependencies: *tidyverse, gsubfn, dummies, data.table, gbm*


## Installation

First, clone the source files from GitHub and install the required packages. 
```
git clone https://github.com/phansol/gremlin.git
```
Then, download `gremlin.fit.rds` and `data.tar.xz` (should be decompressed before use) from [here](ftp_server_address) to `gremlin/`

## Demo
To demo GREMLIN, run the following command (expected run time: 3 minutes). Expected outputs are in `demo/expected_output`
```
gremlin -v demo/input/somatic.svs.callset.sort -t demo/input/tumor.bam -n demo/input/normal.bam -r demo/input/reference.fa -i demo
```

## Usage
```
gremlin [-h] [-v CALL_SET] [-n NORMAL_BAM] [-t TUMOR_BAM] [-r REFERENCE_FASTA]
        [-i SAMPLE_ID] [-o OUTPUT_DIRECTORY] [-g REFERENCE_VERSION]
        [-c TUMOR_CELL_FRACTION] [-p TUMOR_PLOIDY] [-w WGD_STATUS] [-y TUMOR_TISSUE]
```
#### Arguments:
##### &nbsp;&nbsp;&nbsp; (Required)
* ``-v`` Structural variation call set (tab-delimited; with a header CHR1 POS1 CHR2 POS2 SVTYPE CT <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(SVTYPE: DEL|DUP|INV|TRA, CT: 3to3|3to5|5to3|5to5)
* ``-n`` Normal bam (or cram)
* ``-t`` Tumor bam (or cram)
* ``-r`` Reference fasta (index should be [given_fasta].fai)

##### &nbsp;&nbsp;&nbsp; (Optional)
* ``-i`` Sample ID [default: basename of input call set]
* ``-o`` Output directory [default: directory of input call set]
* ``-g`` Reference version (19|38) [default: 19]
* ``-c`` Tumor cell fraction [default: 0.5]
* ``-p`` Tumor genome ploidy [default: 2]
* ``-w`` Whole-genome duplication status (wgd|no_wgd) [default: no_wgd]
* ``-y`` Tumor tissue (Biliary|Bladder|Bone_SoftTissue|Breast|Cervix|CNS|Colon_Rectum|Esophagus|Head_Neck|<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Hematologic|Kideny|Liver|Lung|Ovary|Pancreas|Prostate|Skin|Stomach|Thyroid|Uterus) [default: Biliary] 

#### Output:
* ``[OUTPUT_DIRECTORY]/[SAMPLE_ID].gremlin.feature.dummies.pon.score``: scored SV call set

* ``[OUTPUT_DIRECTORY]/[SAMPLE_ID].gremlin.somatic.svs.tier1.tsv``: tier1 SV calls predicted to be true somatic mutations by GREMLIN

* ``[OUTPUT_DIRECTORY]/[SAMPLE_ID].gremlin.somatic.svs.tier2.tsv``: tier2 SV calls refined with more lenient filtering threshold than tier1

## Best practice
|Step|Description|
|:--:|--|
|1|*(Optional)* [Quality control of input sequences](#quality-control-of-input-sequences)|
|2|[Preprocessing of input SV call sets](#preprocessing-of-input-sv-call-sets)|
|3|[Applying GREMLIN](#usage)|
|4|*(Optional)* [Adjusting classification threshold](#adjusting-classification-threshold)<br>*(Optional)* [Retraining GREMLIN](#retraining-gremlin)<br>*(Optional)* [Additional filtering using normal panels of your cohort](#additional-filtering-using-normal-panels-of-your-cohort)|

## Quality control of input sequences
Short inversion artifacts and artificial fluctuations in sequencing coverage are major sources of false-positive SV calls, commonly seen in whole-genome sequences of low-quality genomic DNA. Thus, we recommend checking your sequencing data as follows before applying GREMLIN. 

Before running the following commands, install required packages using `Rscript codes/requirements.qc.R`

### Flag for short inversion artifacts
The following command will estimate the fraction of short inversions among total read pairs using [samtools](http://www.htslib.org/). If your data has an exceptionally high fraction of short inversions, you will get a fail flag, and the refined list (GREMLIN’s output) may include many short inversion errors.
```
Usage: sh codes/1_quality_check/samtools_short_inv.sh [TUMOR_BAM/CRAM] [THREADS]

Output: [TUMOR_BAM/CRAM].shinv.pass or [TUMOR_BAM/CRAM].shinv.fail
```

### Flag for variable sequencing coverage
The artificial fluctuations in the sequencing coverage can be estimated as follows. The number of aligned reads across the genome will be inferred by indexcov which can be installed from [here](https://github.com/brentp/goleft/tree/master/indexcov). Then, the depth ratio between tumor and normal sequences will be fitted into step functions to offset the variabilities derived from copy number variations. The overall coverage fluctuation will be measured as the mean squared deviation between the depth ratios and the fitted lines.
```
Usage for bam: sh codes/1_quality_check/indexcov_read_depth.bam.sh [TUMOR_BAM] [NORMAL_BAM] [OUTPUT_DIRECTORY] [REFERENCE_BUILD] [REFERENCE_FASTA_INDEX]
Usage for cram: sh codes/1_quality_check/indexcov_read_depth.cram.sh [TUMOR_CRAM] [NORMAL_CRAM] [OUTPUT_DIRECTORY] [REFERENCE_BUILD] [REFERENCE_FASTA_INDEX]

Output: [TUMOR_BAM/CRAM].depth_ratio.png
        [TUMOR_BAM/CRAM].depth.pass or [TUMOR_BAM/CRAM].depth.fail
```
#### Arguments:
* ``REFERENCE_BUILD`` Reference genome version (19|38)
* ``REFERENCE_FASTA_INDEX`` /path/to/reference.fasta.fai


## Preprocessing of input SV call sets
If you called SVs using DELLY, SvABA, BRASS, or dRanger, run the following command.
```
Usage: Rscript codes/2_preprocessing/vcf_formatting.R [VCF] [CALLER] [REFERENCE_FASTA_INDEX] [OUTPUT_DIRECTORY]

Output: [OUTPUT_DIRECTORY]/[VCF].sort
```
#### Arguments:
* ``VCF`` VCF file obtained from DELLY, SvABA, BRASS, or dRanger
* ``CALLER`` SV caller used to get the call set (DELLY|SvABA|BRASS|dRanger)

Otherwise, transform your SV call set into the following tab-separated format. 
```
(e.g.,)	CHR1	POS1	CHR2	POS2	SVTYPE	CT
	1	6435211	1	6475356	DEL	3to5
	3	2762555	3	3546843	DUP	5to3
	8	154337	8	165439	INV	3to3
	10	1293850	X	3287959	TRA	5to3	
```

## Adjusting classification threshold
You can adjust the filtering threshold (default is 0.89 for tier1 and 0.57 for tier2).
```
Usage: Rscript codes/5_postprocessing/optional_adjusting_classification_threshold.R [OUTPUT] [THRESHOLD]

Output: [SAMPLE_ID].gremlin.somatic.svs.cutoff_[THRESHOLD].tsv
```
#### Arguments:
* ``OUTPUT`` GREMLIN's output *([SAMPLE_ID].gremlin.feature.dummies.pon.score)*
* ``THRESHOLD`` Classification threshold (between 0 and 1)

## Retraining GREMLIN
Install R dependencies using `Rscript codes/requirements.rt.R`

### 1. Retraining with your data
```
Usage: Rscript codes/5_postprocessing/optional_retraining.R [NEW_DATASET] [OUTPUT_DIRECTORY] [SUFFIX] [PERCENT]

Output: [OUTPUT_DIRECTORY]/training_set_[SUFFIX].tsv
	[OUTPUT_DIRECTORY]/gremlin_retrained_[SUFFIX].fit.rds
```
#### Arguments:
* ``NEW_DATASET`` Feature-annotated call set to include for retraining GREMLIN <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Same format as *\*.gremlin.feature.dummies.pon.score* with an additional column "true_label" = T or F
* ``PERCENT`` *(Optional)* The percent of our training set to be included in model retraining (between 0 and 100) <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;160 training samples will be used in default
	       
### 2. Applying the retrained model to your data
```
Usage: Rscript codes/5_postprocessing/optional_apply_retrained_gremlin.R [OUTPUT] [re_trained_gbm.fit.rds] [THRESHOLD]

Output: [OUTPUT].retrained_score
	[SAMPLE_ID].retrained_gremlin.somatic.svs.cutoff_[THRESHOLD].tsv (if threshold is given)
```
#### Arguments:
* ``OUTPUT`` GREMLIN's output *([SAMPLE_ID].gremlin.feature.dummies.pon.score)*
* ``THRESHOLD`` *(Optional)* Classification threshold (between 0 and 1)

## Additional filtering using normal panels of your cohort
```
Usage: Rscript codes/5_postprocessing/optional_cohort_specific_pon_annotation.R [OUTPUT] [PON] [COHORT_ID]

Output: [OUTPUT].pon_[COHORT_ID]
```
#### Arguments:
* ``OUTPUT`` GREMLIN's output (\*.feature.dummies.pon.score)
* ``PON`` Your panel of normal data satisfying the following conditions <br>
	- Column order should be (1) CHR1, (2) POS1, (3) CHR2, (4) POS2, (5) CT, (6) SVTYPE, and (7) SAMPLE_ID <br>
	- Each line should be sorted as CHR1 <= CHR2 and POS1 <= POS2 <br>
	- Tab-separated without column names <br>
	```
	(e.g.,) 2	648899	2	1238794	3to5	DEL	sample_id_1
		 5	6876412	5	7425230	3to5	DEL	sample_id_1
		 7	1215465	22	75212	5to3	TRA	sample_id_2
		 11	4373522	11	4588301	5to5	INV	sample_id_3
	```
* ``COHORT_ID`` Used for the column name of your PoN "normal_panel_\[COHORT_ID]"
