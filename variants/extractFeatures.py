import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import variants
import func
import pandas
import os
import features
from multiprocessing import Pool


print sys.argv
if len(sys.argv) != 3:
    sys.exit('two arguments required: child_id output/dir')
#child_id = '11006.p1'
#caller = 'fb'
child_id = sys.argv[1]
output_directory = sys.argv[2]
min_DP = 7
var_type = 'SNP'
vcf_pat = '/mnt/scratch/asalomatov/data/columbia/vcf/%s_%s-02_%s-01.annotated-norm.vcf.gz'
#vcf_pat = '/mnt/scratch/asalomatov/data/columbia/rerun_fb/%s-norm.vcf.gz'
bam_pat = '/mnt/scratch/asalomatov/data/columbia/bam/Sample.%s.bam'
bai_pat = '/mnt/scratch/asalomatov/data/columbia/bam/Sample.%s.bam.bai'
infile_ped = '/mnt/scratch/asalomatov/data/columbia/pcgc_ped.txt'
known_vars = '/mnt/scratch/asalomatov/data/columbia/pcgc_denovo_snp.tsv'
output_dir = output_directory
output_dir_known = ''
func.runInShell('mkdir -p ' + output_dir)
if known_vars:
    output_dir_known = os.path.join(output_dir, 'known')
    func.runInShell('mkdir -p ' + output_dir_known)


def multi_wrap_readBamReadcount(args):
    return func.readBamReadcount(*args)

### SSC ped
#ped = reload(ped)
#func = reload(func)
#features = reload(features)

myped = ped.Ped(infile_ped)
myped.addVcf(file_pat=vcf_pat)
#myped.addVcf()
myped.ped.dropna(subset=['vcf'], inplace=True)
myped.ped.shape
#myped.addBam()
myped.addBam(file_pat=bam_pat)
myped.ped.dropna(subset=['bam'], inplace=True)
myped.ped.shape
myped.addBai(file_pat=bai_pat)
myped.ped.dropna(subset=['bai'], inplace=True)
myped.ped.shape
myped.ped.head()

f = features.Features(myped, known_vars)
if not f.initTrioFor(child_id):
    sys.stderr.write('\nfailed to initialize trio for ' + child_id)
else:
    sys.stdout.write('\ninitialized trio for ' + child_id)
    f.extractFeatures()
    pool = Pool(3)
    fam_features = pool.map(multi_wrap_readBamReadcount,
                            [(f.sample_features,),
                             (f.father_features,),
                             (f.mother_features,)])
    f.removeTmpDir()
    fam_features[0].columns = ['CHROM', 'POS'] + func.addSuffix(fam_features[0].columns[2:], '_offspring')
    fam_features[1].columns = ['CHROM', 'POS'] + func.addSuffix(fam_features[1].columns[2:], '_father')
    fam_features[2].columns = ['CHROM', 'POS'] + func.addSuffix(fam_features[2].columns[2:], '_mother')
    fam_f = fam_features[0]
    fam_f = fam_f.merge(fam_features[1], how='inner', on=['CHROM', 'POS'])
    fam_f = fam_f.merge(fam_features[2], how='inner', on=['CHROM', 'POS'])
    fam_f['ind_id'] = child_id
    # mark verified if information is available
    if known_vars:
        f.verified_variants.columns = ['ind_id', 'CHROM', 'POS', 'ref', 'alt', 'status',  'descr']
        f.verified_variants['POS'] = f.verified_variants['POS'].astype(str)
        fam_f = fam_f.merge(f.verified_variants[['ind_id', 'CHROM', 'POS', 'status', 'descr']],
                            how='left', on=['ind_id', 'CHROM', 'POS'])
        fam_f.ix[fam_f.status.isnull(), 'status'] = 'NA'
    clmns_to_drop = [# 'REF_offspring',
                     'REF_father',
                     'REF_mother',
                     'DP_offspring',
                     'DP_mother',
                     'DP_father',
                     'REF_avg_se_mapping_quality_offspring',
                     'REF_avg_se_mapping_quality_father',
                     'REF_avg_se_mapping_quality_mother',
                     'ALT_avg_se_mapping_quality_offspring',
                     'ALT_avg_se_mapping_quality_father',
                     'ALT_avg_se_mapping_quality_mother',
                     'INDEL_avg_se_mapping_quality_offspring',
                     'INDEL_avg_se_mapping_quality_father',
                     'INDEL_avg_se_mapping_quality_mother',
                     # 'ALT_base_offspring',
                     'ALT_base_father',
                     'ALT_base_mother',
                     'INDEL_base_offspring',
                     'INDEL_base_father',
                     'INDEL_base_mother']
    fam_f.drop(clmns_to_drop, axis=1, inplace=True)
    fam_f['DP_offspring'] = fam_f['REF_count_offspring'] + fam_f['ALT_count_offspring']
    fam_f['DP_father'] = fam_f['REF_count_father'] + fam_f['ALT_count_father']
    fam_f['DP_mother'] = fam_f['REF_count_mother'] + fam_f['ALT_count_mother']
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring)
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring != 0)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father != 0)
    # filter out DP < 8
    c1 = fam_f['DP_offspring'] > min_DP
    c2 = fam_f['DP_father'] > min_DP
    c3 = fam_f['DP_mother'] > min_DP
    fam_f = fam_f[ c1 & c2 & c3 ]
    fam_f.DP_offspring.value_counts()
    fam_f.columns
    fam_f.dtypes
    print fam_f.shape
    fam_f.to_csv(os.path.join(output_dir, child_id), sep='\t', index=False)
    x = fam_f[~fam_f.status.isin(['NA'])]
    if not x.empty:
        x.to_csv(os.path.join(output_dir_known, child_id), sep='\t', index=False)






