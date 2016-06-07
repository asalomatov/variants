input_tsv=$1
script_dir=$2
child_id=$3
work_dir="$(dirname $input_tsv)"
base_name=$(basename "$input_tsv" | cut -d. -f1)
echo $base_name
#cd $work_dir

cp ${script_dir}/header_extra.txt ${work_dir}/${child_id}.vcf
cat $input_tsv | sort -V -k1,1 -k2,2 >> ${work_dir}/${child_id}.vcf
make -f ${script_dir}/annSnpEff.mk PREFIX=${child_id} SUFFIX=.vcf INDIR=$work_dir OUTDIR=$work_dir
cat ${work_dir}/${child_id}-ann.vcf | \
    /bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/share/scripts/vcfEffOnePerLine.pl > \
    ${work_dir}/${child_id}-ann-onePline.vcf
${script_dir}/extractFields.sh ${work_dir}/${child_id}-ann-onePline.vcf > ${work_dir}/${child_id}-ann-onePline.tsv
