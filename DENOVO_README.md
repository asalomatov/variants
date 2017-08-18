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

Install **_variants_** from github

```
pip install git+git://github.com/asalomatov/variants.git
```

Install prerequisite 

```
pip install tensorflow
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
ind_id,CHROM,POS,REF,ALT,pred_prob,DP_of,DP_fa,DP_mo,var_id,var_id_a,alleles_of,alleles_fa,alleles_mo,num_alt_of,num_alt_fa,num_alt_mo,num_alt_all,inherit_fa,inherit_mo,inherit_prnts
trio001.p1,1,12073513,C,T,0.834789872169,10,28,26,trio001.p1_1_12073513,trio001.p1_1_12073513_C_T,C_8_T_2,C_28_none,C_26_none,1,0,0,1,0,0,0
trio001.p1,1,12939426,T,C,2.73509740509e-05,49,135,162,trio001.p1_1_12939426,trio001.p1_1_12939426_T_C,T_40_C_9,T_107_C_27_G_1,T_162_none,1,2,0,3,1,0,1
trio001.p1,1,12939440,C,G,0.000129409629153,49,142,156,trio001.p1_1_12939440,trio001.p1_1_12939440_C_G,C_38_G_11,C_116_A_1_G_25,C_153_A_3,1,2,1,4,1,0,1
trio001.p1,1,12939476,G,C,6.94080290486e-06,57,139,172,trio001.p1_1_12939476,trio001.p1_1_12939476_G_C,G_37_C_20,G_135_C_4,G_120_C_52,1,1,1,3,1,1,2
```
