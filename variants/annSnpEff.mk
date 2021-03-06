### 
SHELL = /bin/bash
USR = $(shell whoami)
INCLMK = ~/projects/pipeline/ppln/include.mk
include $(INCLMK)
### may override on cl
PREFIX = 1
SUFFIX = 
INDIR = .
OUTDIR = .
LOGDIR = $(OUTDIR)
TMPDIR = /tmp/$(USR)
#SNPSIFTJAR = /mnt/xfs1/bioinfoCentos7/software/installs/bcbio_nextgen/150617/Cellar/snpeff/4.1g/libexec/SnpSift.jar
#SNPEFFJAR = /mnt/xfs1/bioinfoCentos7/software/installs/bcbio_nextgen/150617/Cellar/snpeff/4.1g/libexec/snpEff.jar
VEP = /mnt/xfs1/bioinfoCentos7/software/builds/perlbrew/cellar/perls/perl-5.22.1/bin/perl /mnt/xfs1/bioinfoCentos7/software/builds/ensembl/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl --cache --dir=/mnt/ceph/users/carriero/VEPcache/.vep --everything --assembly $(VEPGENBUILD) --offline --force_overwrite -o STDOUT --fork 4 --fasta $(GENOMEREF) --no_stats

###
inFile = $(wildcard $(INDIR)/$(PREFIX)*$(SUFFIX))
$(info $(inFile))
outFile = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX),%-ann$(SUFFIX),$(notdir $(inFile))))
sumFile = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX),%-ann$(SUFFIX).summary.html,$(notdir $(inFile))))
$(info $(outFile))

all: $(outFile)

$(outFile): $(inFile)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR) annotate -c $(SNPEFFCONF) $(DBSNP) -id -noLog $< | \
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR) varType  -c $(SNPEFFCONF) -noLog - | \
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR) annotate  -c $(SNPEFFCONF) -noLog -tabix $(DBSPIDEX) - | \
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR)  dbnsfp - -a -c $(SNPEFFCONF) -db $(DBNSFP) \
	-f rs_dbSNP147,aapos,aaref,aaalt,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,PROVEAN_score,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,M-CAP_score,M-CAP_rankscore,M-CAP_pred,1000Gp3_AF,ExAC_AF,ALSPAC_AC,ALSPAC_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF | sed "s/dbNSFP_GERP++/dbNSFP_GERP/g" | sed "s/dbNSFP_M-CAP/dbNSFP_M_CAP/g" | \
	$(JAVA) -Xmx5G -jar $(SNPEFFJAR) ann -nextProt -noLog -noStats -c $(SNPEFFCONF) $(SNPEFFGENOME) -v -lof > $@
	cat $@ | $(SNPEFFDIR)/scripts/vcfEffOnePerLine.pl > $@.onePline.vcf
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR) extractFields -s "__" -e "ZZZ" -noLog $@.onePline.vcf  CHROM POS ID REF ALT QUAL VARTYPE \
    "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" \
    "ANN[*].FEATURE" "ANN[*].FEATUREID" \
    "ANN[*].BIOTYPE" "Ann[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \
    "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" \
    "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS" \
    "LOF[*].GENE" "LOF[*].GENEID" "LOF[*].NUMTR" "LOF[*].PERC" \
    "NMD[*].GENE" "NMD[*].GENEID" "NMD[*].NUMTR" "NMD[*].PERC" \
    "ind_id" "pred_labels" "pred_prob" "status" "ref_DP" "alt_DP" \
    "DP" "DP_of" "DP_fa" "DP_mo" \
    "alleles_of" "alleles_fa" "alleles_mo" \
    "num_alt_of" "num_alt_fa" "num_alt_mo" "num_alt_all" \
    "inherit_fa" "inherit_mo" "inherit_prnts" \
    "var_id" "var_id_a" \
    "dbNSFP_PROVEAN_score" "dbNSFP_CADD_raw_rankscore" \
    "dbNSFP_GERP_RS"  "dbNSFP_genename"  \
    "dbNSFP_GERP_NR" "dbNSFP_Ensembl_transcriptid"  \
    "dbNSFP_Polyphen2_HVAR_score"  "dbNSFP_MutationAssessor_score" \
    "dbNSFP_Ensembl_proteinid" "dbNSFP_rs_dbSNP147" \
    "dbNSFP_aapos" \
    "dbNSFP_aaref" \
    "dbNSFP_aaalt" \
    "dbNSFP_Uniprot_acc_Polyphen2" \
    "dbNSFP_Uniprot_id_Polyphen2" \
    "dbNSFP_Uniprot_aapos_Polyphen2" \
    "dbNSFP_Ensembl_geneid" "dbNSFP_CADD_raw" \
    "dbNSFP_GERP_RS_rankscore" \
    "dbNSFP_SIFT_score" \
    "dbNSFP_1000Gp3_AF" \
    "dbNSFP_ExAC_AF" "dbNSFP_ALSPAC_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF" \
    "dbNSFP_Polyphen2_HDIV_score" \
    "dbNSFP_PROVEAN_pred" \
    "dbNSFP_CADD_phred" \
    "dbNSFP_Polyphen2_HDIV_pred" \
    "dbNSFP_MutationTaster_score" \
    "dbNSFP_MutationTaster_pred" \
    "dbNSFP_MutationAssessor_pred" \
    "dbNSFP_Polyphen2_HVAR_pred" \
    "dbNSFP_SIFT_pred" \
    "dbNSFP_MetaSVM_score" \
    "dbNSFP_MetaSVM_rankscore" \
    "dbNSFP_MetaSVM_pred" \
    "dbNSFP_MetaLR_score" \
    "dbNSFP_MetaLR_rankscore" \
    "dbNSFP_MetaLR_pred" \
    "dbNSFP_M_CAP_score" "dbNSFP_M_CAP_rankscore" "dbNSFP_M_CAP_pred" \
    "spidex_dpsi_max_tissue" \
    "spidex_dpsi_zscore" \
    "spidex_gene" \
    "spidex_strand" \
    "spidex_transcript" \
    "spidex_exon_number" \
    "spidex_location" \
    "spidex_cds_type" \
    "spidex_ss_dist"	> $@.onePline.tsv


#	$(VEP) --vcf --output_file STDOUT --input_file $< | \
#	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR) annotate -id -c $(SNPEFFCONF) -noLog -dbsnp $< | \
#       $(JAVA) -Xmx5G -jar $(SNPSIFTJAR) varType -c $(SNPEFFCONF) -noLog - | \
#	$(JAVA) -Xmx5G -jar $(SNPEFFJAR) ann -noLog -noStats -c $(SNPEFFCONF) $(SNPEFFGENOME) -v -lof  | \


##-csvStats -s $(sumFile)
# '
# chr
# pos(1-coor)
# ref
# alt
# aaref
# aaalt
# rs_dbSNP147
# hg19_chr
# hg19_pos(1-coor)
# hg18_chr
# hg18_pos(1-coor)
# genename
# cds_strand
# refcodon
# codonpos
# codon_degeneracy
# Ancestral_allele
# AltaiNeandertal
# Denisova
# Ensembl_geneid
# Ensembl_transcriptid
# Ensembl_proteinid
# aapos
# SIFT_score
# SIFT_converted_rankscore
# SIFT_pred
# Uniprot_acc_Polyphen2
# Uniprot_id_Polyphen2
# Uniprot_aapos_Polyphen2
# Polyphen2_HDIV_score
# Polyphen2_HDIV_rankscore
# Polyphen2_HDIV_pred
# Polyphen2_HVAR_score
# Polyphen2_HVAR_rankscore
# Polyphen2_HVAR_pred
# LRT_score
# LRT_converted_rankscore
# LRT_pred
# LRT_Omega
# MutationTaster_score
# MutationTaster_converted_rankscore
# MutationTaster_pred
# MutationTaster_model
# MutationTaster_AAE
# MutationAssessor_UniprotID
# MutationAssessor_variant
# MutationAssessor_score
# MutationAssessor_score_rankscore
# MutationAssessor_pred
# FATHMM_score
# FATHMM_converted_rankscore
# FATHMM_pred
# PROVEAN_score
# PROVEAN_converted_rankscore
# PROVEAN_pred
# Transcript_id_VEST3
# Transcript_var_VEST3
# VEST3_score
# VEST3_rankscore
# MetaSVM_score
# MetaSVM_rankscore
# MetaSVM_pred
# MetaLR_score
# MetaLR_rankscore
# MetaLR_pred
# Reliability_index
# M-CAP_score
# M-CAP_rankscore
# M-CAP_pred
# CADD_raw
# CADD_raw_rankscore
# CADD_phred
# DANN_score
# DANN_rankscore
# fathmm-MKL_coding_score
# fathmm-MKL_coding_rankscore
# fathmm-MKL_coding_pred
# fathmm-MKL_coding_group
# Eigen_coding_or_noncoding
# Eigen-raw
# Eigen-phred
# Eigen-PC-raw
# Eigen-PC-phred
# Eigen-PC-raw_rankscore
# GenoCanyon_score
# GenoCanyon_score_rankscore
# integrated_fitCons_score
# integrated_fitCons_score_rankscore
# integrated_confidence_value
# GM12878_fitCons_score
# GM12878_fitCons_score_rankscore
# GM12878_confidence_value
# H1-hESC_fitCons_score
# H1-hESC_fitCons_score_rankscore
# H1-hESC_confidence_value
# HUVEC_fitCons_score
# HUVEC_fitCons_score_rankscore
# HUVEC_confidence_value
# GERP++_NR
# GERP++_RS
# GERP++_RS_rankscore
# phyloP100way_vertebrate
# phyloP100way_vertebrate_rankscore
# phyloP20way_mammalian
# phyloP20way_mammalian_rankscore
# phastCons100way_vertebrate
# phastCons100way_vertebrate_rankscore
# phastCons20way_mammalian
# phastCons20way_mammalian_rankscore
# SiPhy_29way_pi
# SiPhy_29way_logOdds
# SiPhy_29way_logOdds_rankscore
# 1000Gp3_AC
# 1000Gp3_AF
# 1000Gp3_AFR_AC
# 1000Gp3_AFR_AF
# 1000Gp3_EUR_AC
# 1000Gp3_EUR_AF
# 1000Gp3_AMR_AC
# 1000Gp3_AMR_AF
# 1000Gp3_EAS_AC
# 1000Gp3_EAS_AF
# 1000Gp3_SAS_AC
# 1000Gp3_SAS_AF
# TWINSUK_AC
# TWINSUK_AF
# ALSPAC_AC
# ALSPAC_AF
# ESP6500_AA_AC
# ESP6500_AA_AF
# ESP6500_EA_AC
# ESP6500_EA_AF
# ExAC_AC
# ExAC_AF
# ExAC_Adj_AC
# ExAC_Adj_AF
# ExAC_AFR_AC
# ExAC_AFR_AF
# ExAC_AMR_AC
# ExAC_AMR_AF
# ExAC_EAS_AC
# ExAC_EAS_AF
# ExAC_FIN_AC
# ExAC_FIN_AF
# ExAC_NFE_AC
# ExAC_NFE_AF
# ExAC_SAS_AC
# ExAC_SAS_AF
# ExAC_nonTCGA_AC
# ExAC_nonTCGA_AF
# ExAC_nonTCGA_Adj_AC
# ExAC_nonTCGA_Adj_AF
# ExAC_nonTCGA_AFR_AC
# ExAC_nonTCGA_AFR_AF
# ExAC_nonTCGA_AMR_AC
# ExAC_nonTCGA_AMR_AF
# ExAC_nonTCGA_EAS_AC
# ExAC_nonTCGA_EAS_AF
# ExAC_nonTCGA_FIN_AC
# ExAC_nonTCGA_FIN_AF
# ExAC_nonTCGA_NFE_AC
# ExAC_nonTCGA_NFE_AF
# ExAC_nonTCGA_SAS_AC
# ExAC_nonTCGA_SAS_AF
# ExAC_nonpsych_AC
# ExAC_nonpsych_AF
# ExAC_nonpsych_Adj_AC
# ExAC_nonpsych_Adj_AF
# ExAC_nonpsych_AFR_AC
# ExAC_nonpsych_AFR_AF
# ExAC_nonpsych_AMR_AC
# ExAC_nonpsych_AMR_AF
# ExAC_nonpsych_EAS_AC
# ExAC_nonpsych_EAS_AF
# ExAC_nonpsych_FIN_AC
# ExAC_nonpsych_FIN_AF
# ExAC_nonpsych_NFE_AC
# ExAC_nonpsych_NFE_AF
# ExAC_nonpsych_SAS_AC
# ExAC_nonpsych_SAS_AF
# clinvar_rs
# clinvar_clnsig
# clinvar_trait
# clinvar_golden_stars
# Interpro_domain
# GTEx_V6_gene
# GTEx_V6_tissue
# '
