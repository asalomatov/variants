#!/bin/env python
from __future__ import print_function
import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import variants
import func
import pandas
import os
import features
import yaml
from multiprocessing import Pool

print(sys.argv)
if len(sys.argv) != 3:
    sys.exit('two arguments required: child_id path/to/config.yml')

child_id = sys.argv[1]
config_file = sys.argv[2]

# get parameters from yaml config file
with open(config_file, 'r') as f:
    cfg = yaml.safe_load(f)

min_DP = cfg['min_DP']
var_type = cfg['variant_type']
vcf_pat = cfg['vcf_pattern']
bam_pat = cfg['bam_pattern']
bai_pat = bam_pat + '.bai'
ped_file = cfg['ped_file']
ped_file_extended = cfg['ped_file_extended']
bam_readcount = cfg['bam_readcount']
genome_ref = cfg['genome_ref']
known_vars = cfg['known_variants']
output_dir = cfg['output_directory']
#output_dir_known = ''

# create output dirs
func.runInShell('mkdir -p ' + output_dir)
if known_vars:
    output_dir_known = os.path.join(output_dir, 'known')
    func.runInShell('mkdir -p ' + output_dir_known)


# wrap a funtion for use with multiprocessing
def multi_wrap_readBamReadcount(args):
    return func.readBamReadcount(*args)

# populate ped DF
myped = ped.Ped(ped_file)
myped.addVcf(file_pat=vcf_pat)
myped.ped.dropna(subset=['vcf'], inplace=True)
myped.addBam(file_pat=bam_pat)
myped.ped.dropna(subset=['bam'], inplace=True)
myped.addBai(file_pat=bai_pat)
myped.ped.dropna(subset=['bai'], inplace=True)

f = features.Features(myped, known_vars)
# trio has to be complete with no files missing
if not f.initTrioFor(child_id):
    sys.stderr.write('\nfailed to initialize trio for ' + child_id)
else:
    sys.stdout.write('\ninitialized trio for ' + child_id)
    sys.stdout.write('\n')
    f.extractFeatures(genome_ref, bam_readcount, var_type)
    pool = Pool(3)
    fam_features = pool.map(multi_wrap_readBamReadcount,
                            [(f.sample_features, var_type),
                             (f.father_features, var_type),
                             (f.mother_features, var_type)])
    f.removeTmpDir()
    fam_features[0].columns = ['CHROM', 'POS'] + func.addSuffix(fam_features[0].columns[2:], '_offspring')
    fam_features[1].columns = ['CHROM', 'POS'] + func.addSuffix(fam_features[1].columns[2:], '_father')
    fam_features[2].columns = ['CHROM', 'POS'] + func.addSuffix(fam_features[2].columns[2:], '_mother')
    fam_f = fam_features[0]
    fam_f = fam_f.merge(fam_features[1], how='inner', on=['CHROM', 'POS'])
    fam_f = fam_f.merge(fam_features[2], how='inner', on=['CHROM', 'POS'])
    fam_f = fam_f[fam_f.ALT_count_offspring.astype(int) != 0]  # remove not-variant loci
    fam_f['ind_id'] = child_id
    # mark verified if information is available
    if known_vars:
        f.verified_variants.columns = ['ind_id', 'CHROM', 'POS', 'ref', 'alt', 'status',  'descr']
        if var_type.lower() == 'indel':
            # deletions should have POS + 1
            c1 = f.verified_variants.ref.apply(len) > 1  # deletion 1
            c2 = f.verified_variants.alt.apply(len) == 1  # deletion 2
            f.verified_variants.ix[:, 'POS'] = f.verified_variants.POS.astype(int)
            f.verified_variants.ix[c1 & c2, 'POS'] = f.verified_variants.POS[c1 & c2] + 1
            f.verified_variants['POS'] = f.verified_variants['POS'].astype(str)
        fam_f = fam_f.merge(f.verified_variants[['ind_id', 'CHROM', 'POS', 'status', 'descr']],
                            how='left', on=['ind_id', 'CHROM', 'POS'])
        fam_f.ix[fam_f.status.isnull(), 'status'] = 'NA'
    clmns_to_drop = ['REF_avg_se_mapping_quality_offspring',
                     'REF_avg_se_mapping_quality_father',
                     'REF_avg_se_mapping_quality_mother',
                     'ALT_avg_se_mapping_quality_offspring',
                     'ALT_avg_se_mapping_quality_father',
                     'ALT_avg_se_mapping_quality_mother',
                     'INDEL_avg_se_mapping_quality_offspring',
                     'INDEL_avg_se_mapping_quality_father',
                     'INDEL_avg_se_mapping_quality_mother',
                     # 'ALT_base_offspring',
                     # 'ALT_base_father',
                     # 'ALT_base_mother',
                     'INDEL_base_offspring',
                     'INDEL_base_father',
                     'INDEL_base_mother']
    # fam_f.drop(clmns_to_drop, axis=1, inplace=True)
    fam_f['DP_offspring'] = fam_f['REF_count_offspring'] + fam_f['ALT_count_offspring']
    fam_f['DP_father'] = fam_f['REF_count_father'] + fam_f['ALT_count_father']
    fam_f['DP_mother'] = fam_f['REF_count_mother'] + fam_f['ALT_count_mother']
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring)
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring != 0)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father != 0)
    # filter out DP < 8
    c1 = fam_f['DP_offspring'] > min_DP - 1
    c2 = fam_f['DP_father'] > min_DP - 1
    c3 = fam_f['DP_mother'] > min_DP - 1
    fam_f = fam_f[c1 & c2 & c3]
    fam_f.to_csv(os.path.join(output_dir, child_id), sep='\t', index=False)
    if known_vars:
        x = fam_f[~fam_f.status.isin(['NA'])]
        if not x.empty:
            x.to_csv(os.path.join(output_dir_known, child_id), sep='\t', index=False)






