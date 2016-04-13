import sys
#sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import variants
import func
import pandas
import os
import features
from multiprocessing import Pool


print sys.argv
if len(sys.argv) != 4:
    sys.exit('three arguments required: child_id caller output/dir')
#child_id = '11006.p1'
#caller = 'fb'
child_id = sys.argv[1]
caller = sys.argv[2]
output_directory = sys.argv[3]
min_DP = 7
var_type = 'SNP'
vcf_pat = ''
if caller == 'hc':
    vcf_pat = '/mnt/scratch/asalomatov/data/SSC/wes/vcf/hc/%s.family.vqsr.sorted-norm.vcf'
elif caller == 'fb':
    vcf_pat = '/mnt/scratch/asalomatov/data/SSC/wes/vcf/fb/%s.family.freebayes.sorted-norm.vcf'
else:
    sys.exit('unknown caller, exiting... ')

output_dir = os.path.join(output_directory, caller)
all_dir = output_dir + '/all_' + var_type
known_dir = output_dir + '/known_' + var_type
func.runInShell('mkdir -p ' + all_dir)
func.runInShell('mkdir -p ' + known_dir)

def multi_wrap_readBamReadcount(args):
    return func.readBamReadcount(*args)



### SSC ped
#ped = reload(ped)
#func = reload(func)
#features = reload(features)

infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
myped = ped.Ped(infile_ped, ['collection'])
myped.addVcf(file_pat=vcf_pat)
#myped.addVcf()
myped.ped.dropna(subset=['vcf'], inplace=True)
myped.ped.shape
#myped.addBam()
myped.addBam(file_pat='/mnt/scratch/asalomatov/data/SSC/wes/bam/%s.realigned.recal.bam')
myped.ped.dropna(subset=['bam'], inplace=True)
myped.ped.shape
myped.addBai(file_pat='/mnt/scratch/asalomatov/data/SSC/wes/bam/%s.realigned.recal.bam.bai')
myped.ped.dropna(subset=['bai'], inplace=True)
myped.ped.shape
myped.ped.head()


#func = reload(func)
#ped = reload(ped)
#features = reload(features)
#variants = reload(variants)

f = features.Features(myped, '/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_denovo_clean_snp.tsv')
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


# add status, and extra fields then save to file
    f.verified_variants.columns = ['ind_id', 'CHROM', 'POS', 'ref', 'alt', 'status',  'descr', 'vartype']
    f.verified_variants['POS'] = f.verified_variants['POS'].astype(str)
    fam_f = fam_f.merge(f.verified_variants[['ind_id', 'CHROM', 'POS', 'status',  'descr', 'vartype']], 
            how='left', on=['ind_id', 'CHROM', 'POS'])
    fam_f.ix[fam_f.status.isnull(), 'status'] = 'NA'
    fam_f.shape
    fam_f.dtypes
    fam_f.status.value_counts()
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
                     'INDEL_base_mother',
                     'vartype']
#clmns_to_keep = [i for i in fam_f if i not in clmns_to_drop]
    fam_f.drop(clmns_to_drop, axis=1, inplace=True)
    fam_f['DP_offspring'] =  fam_f['REF_count_offspring'] + fam_f['ALT_count_offspring']
    fam_f['DP_father'] =  fam_f['REF_count_father'] + fam_f['ALT_count_father']
    fam_f['DP_mother'] =  fam_f['REF_count_mother'] + fam_f['ALT_count_mother']
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring)
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring != 0)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father != 0)
## filter out DP < 8
    c1 = fam_f['DP_offspring'] > min_DP
    c2 = fam_f['DP_father'] > min_DP
    c3 = fam_f['DP_mother'] > min_DP
    fam_f = fam_f[ c1 & c2 & c3 ]
    fam_f.DP_offspring.value_counts()
    fam_f.columns
    fam_f.dtypes
    print fam_f.shape
    fam_f.to_csv(os.path.join(all_dir, child_id), sep='\t', index=False)
    x = fam_f[~fam_f.status.isin(['NA'])]
    if not x.empty:
        x.to_csv(os.path.join(known_dir, child_id), sep='\t', index=False)






