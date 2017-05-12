#!/bin/bash

java -jar /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar extractFields -s "__" -e "ZZZ"  -noLog $1 \
    CHROM POS ID REF ALT QUAL VARTYPE \
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
    "inherit_fa" "inherit_mo" "inherit_all" \
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
    "spidex_ss_dist"
