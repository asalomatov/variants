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
o0 = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX),%-norm.vcf,$(notdir $(inFiles))))
$(info $(o0))

all: $(o0)

$(OUTDIR)/%-norm.vcf: $(INDIR)/%$(SUFFIX)
	mkdir -p $(LOGDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	vt normalize -r $(REF) -o $@ $<
	bgzip -f $@
	tabix -p vcf $@.gz
