
input_tsv=$1
script_dir=$2
prefix=$3
target_bed=$4
incl_make=$5
vep_refseq=$6
work_dir="$(dirname $input_tsv)"
# base_name=$(basename "$input_tsv" | cut -d. -f1)
# echo $base_name
# cd $work_dir

cp ${script_dir}/header_extra.txt ${work_dir}/tmp-${prefix}.vcf
cat $input_tsv | sort -V -k1,1 -k2,2 >> ${work_dir}/tmp-${prefix}.vcf
vcfintersect -b $target_bed ${work_dir}/tmp-${prefix}.vcf > ${work_dir}/${prefix}.vcf
echo 'running VEP'
make -f ${script_dir}/annVEP.mk INCLMK=$incl_make VEPREFSEQ=$vep_refseq PREFIX=${prefix} SUFFIX=.vcf INDIR=$work_dir OUTDIR=$work_dir
echo 'running snpEff'
make -f ${script_dir}/annSnpEff.mk INCLMK=$incl_make PREFIX=${prefix} SUFFIX=.vcf INDIR=$work_dir OUTDIR=$work_dir
cat ${work_dir}/${prefix}-ann.vcf | \
    /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/share/scripts/vcfEffOnePerLine.pl > \
    ${work_dir}/${prefix}-ann-onePline.vcf
${script_dir}/extractFields.sh ${work_dir}/${prefix}-ann-onePline.vcf > ${work_dir}/${prefix}-ann-onePline.tsv
