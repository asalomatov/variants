### *de novo classifier*

### Overview

**_de novo classifier_** is a tool for identifying *de novo* genomic variants. It is 
similar to [DNMFilter](http://www.ncbi.nlm.nih.gov/pubmed/24618463) in its approach but it is
not limited to using gradient-boosting as a sole predictive model. Please see [slides](https://www.dropbox.com/s/ico6qo6pe0zanqe/denovo_filt_IT_20160520.pptx?dl=0) for details on features, model selection, training set, validation ect.

### Getting Started

#### Inslallation

1. Install [bam-readcount](https://github.com/genome/bam-readcount).

2. Install **_variants_** from github 
```
pip install git+git://github.com/asalomatov/variants.git
```
    
#### Configuration

Put together a yaml config file, see 
[example](https://github.com/asalomatov/variants/blob/master/variants/denovo_classifier_config/cfg.yml).
    
You will need a PED file. A sample file:
trio id | sample id | father's id (0 if missing) | mother's id (0 if missing) | sex 1-male 2-female | phenotype 1-unaffected 2-affected |
------- | --------- | -------------------------- | -------------------------- | ------------------- | --------------------------------- |
trio001 | trio001.fa_641964 |  0 |      0 |      1 |      1 |
trio001 | trio001.mo_641942 |  0 |      0 |      2 |      1 |
trio001 | trio001.p1_641943 |  trio001.fa_641964 |      trio001.mo_641942 |      1 |      2 |

A VCF file, and BAM files are specified via patterns to be used with python string substitution.
A `vcf_pattern` like `/path/to/%s-example.vcf.gz` will be transated to `/path/to/trio1-example.vcf.gz`,
this VCF file must have samples `trio001.fa_641964`, `trio001.mo_641942`, `trio001.p1_641943` in the header.
A `bam_pattern` like `/path/to/%s-example.bam` will be transated to a 

```
/path/to/trio001.fa_641964-example.bam
/path/to/trio001.mo_641942-example.bam
/path/to/trio001.p1_641943-example.bam
```

#### Ready to run

Issue in your terminal
```
call_de_novo.py trio001.p1_641943 /path/to/config/cfg.yml 0.3 
```   

To score all possible *de novo* variants decrease the class probability threashold to 0.

#### Output

     Sample output

```
var_id,pred_prob,test_var_alleles,DP_offspring,DP_father,DP_mother
trio001.p1_641943_10_32975032,0.998314760685,C_13_G_20,33,57,68
trio001.p1_641943_10_49668086,0.505049519446,T_7_A_3,10,11,11
trio001.p1_641943_19_44492568,0.822857483557,A_5_G_7,12,14,23
```
