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
SNPSIFTJAR = /mnt/xfs1/bioinfoCentos7/software/installs/bcbio_nextgen/150617/Cellar/snpeff/4.1g/libexec/SnpSift.jar
SNPEFFJAR = /mnt/xfs1/bioinfoCentos7/software/installs/bcbio_nextgen/150617/Cellar/snpeff/4.1g/libexec/snpEff.jar

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
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR) annotate -id -c $(SNPEFFCONF) -noLog -dbsnp $< | \
        $(JAVA) -Xmx5G -jar $(SNPSIFTJAR) varType -c $(SNPEFFCONF) -noLog - | \
        $(JAVA) -Xmx5G -jar $(SNPSIFTJAR) annotate -c $(SNPEFFCONF) -noLog -tabix $(DBSPIDEX) - | \
	$(JAVA) -Xmx5G -jar $(SNPEFFJAR) ann -noLog -c $(SNPEFFCONF) $(SNPEFFGENOME) -v -lof -csvStats -s $(sumFile) | \
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR)  dbnsfp - -a -c $(SNPEFFCONF) -db $(DBNSFP) \
	-f rs_dbSNP146,aapos,aaref,aaalt,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,PROVEAN_score,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,1000Gp3_AF,ExAC_AF | sed "s/dbNSFP_GERP++/dbNSFP_GERP/g" > $@
