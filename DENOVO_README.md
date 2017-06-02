### *de novo classifier*

### Overview

**_de novo classifier_** is a tool for identifying *de novo* genomic variants. It is
similar to [DNMFilter](http://www.ncbi.nlm.nih.gov/pubmed/24618463) in its approach but it is
not limited to using gradient-boosting as a sole predictive model. Please see [slides](https://www.dropbox.com/s/ico6qo6pe0zanqe/denovo_filt_IT_20160520.pptx?dl=0) for details on features, model selection, training set, validation ect.

### Getting Started

#### Installation

Install [bam-readcount](https://github.com/genome/bam-readcount).

If necessary install [conda](http://conda.pydata.org/miniconda.html)

Create, and activate an environment
```
conda create --name test_env python=2
source activate test_env
```

Install prerequestives

```
conda install matplotlib
conda install numpy
conda install scipy
conda install pandas
conda install scikit-learn
pip install xgboost
```

Install **_variants_** from github

```
pip install git+git://github.com/asalomatov/variants.git
```

#### Configuration

Put together a yaml config file, see
[example](https://github.com/asalomatov/variants/blob/master/variants/denovo_classifier_config/cfg.yml).

You will need a PED file with two additional columns assosiating samples with BAM and VCF files. A sample file(tab delimeted and **headerless**):

trio id | sample id | father's id (0 if missing) | mother's id (0 if missing) | sex 1-male 2-female | phenotype 1-unaffected 2-affected | path to BAM | path to VCF
------- | --------- | -------------------------- | -------------------------- | ------------------- | --------------------------------- | -------------------------- | -------------------------
trio001 | trio001.fa |  0 |      0 |      1 |      1 | ~/bam/trio001.fa.bam  | ~/vcf/trio001.vcf
trio001 | trio001.mo |  0 |      0 |      2 |      1 | ~/bam/trio001.mo.bam  | ~/vcf/trio001.vcf
trio001 | trio001.p1 |  trio001.fa |      trio001.mo |      1 |      2 | ~/bam/trio001.p1.bam  | ~/vcf/trio001.vcf

The VCF file specified in the extended pedigree file must contain genotype information about the trio.
Single sample VCF files are not supported.

#### Ready to run

Issue in your terminal
```
call_de_novo.py trio001.p1 /path/to/config/cfg.yml 0.3
```

To score all possible *de novo* variants decrease the class probability threashold to 0.

#### Output

Sample output

```
var_id,pred_prob,test_var_alleles,DP_offspring,DP_father,DP_mother
trio001.p1_10_32975032,0.998314760685,C_13_G_20,33,57,68
trio001.p1_10_49668086,0.505049519446,T_7_A_3,10,11,11
trio001.p1_19_44492568,0.822857483557,A_5_G_7,12,14,23
```
