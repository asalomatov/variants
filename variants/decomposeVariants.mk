### 
default: all
SHELL = /bin/bash
USR = $(shell whoami)
REF = $(shell echo $$GENOMEREF)
$(info $(REF))
### may override on cl
PREFIX = 1
SUFFIX = .vcf.gz
PROJ = mktest
INDIR = .
OUTDIR = .
LOGDIR = $(OUTDIR)
TMPDIR = /tmp/$(USR)/$(PROJ)
###
inFiles = $(wildcard $(INDIR)/$(FAMCODE)*$(SUFFIX))
$(info $(inFiles))
o0 = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX),%-deco.vcf,$(notdir $(inFiles))))
$(info $(o0))

all: $(o0)

$(OUTDIR)/%-deco.vcf: $(INDIR)/%$(SUFFIX)
	mkdir -p $(LOGDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	zcat $< \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| vt decompose -s - \
	| vt normalize -r $(REF) - \
	| bgzip -c > $@.gz
	tabix -p vcf $@.gz
