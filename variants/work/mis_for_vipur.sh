#!/bin/bash

input_tsv=$1
script_dir=$2
#work_dir="$(dirname $input_tsv)"

cp ${script_dir}/header_extra.txt ${input_tsv}.vcf
cat $input_tsv | sort -V -k1,1 -k2,2 >> ${input_tsv}.vcf
/mnt/xfs1/bioinfoCentos7/software/builds/perlbrew/cellar/perls/perl-5.22.1/bin/perl /mnt/xfs1/bioinfoCentos7/software/builds/ensembl/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl --quiet --cache --dir=/mnt/ceph/users/carriero/VEPcache/.vep --assembly=GRCh37 --port 3337 --force_overwrite --everything --fields=SWISSPROT,TREMBLE --tab -i ${input_tsv}.vcf --output_file=${input_tsv}-vep.tsv
