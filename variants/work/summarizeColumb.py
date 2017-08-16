#!/mnt/xfs1/home/asalomatov/miniconda2/bin/python
import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import train
import numpy
import pandas
import ped
import os



kv_file = '/mnt/scratch/asalomatov/data/columbia/pcgc_denovo_snp.tsv'
kv_vcf_file = '/mnt/scratch/asalomatov/data/columbia/feature_sets/known/all_known.txt'
ped_file = '/mnt/scratch/asalomatov/data/columbia/pcgc_ped.txt'
test_set_pat = '/mnt/scratch/asalomatov/data/columbia/feature_sets/%s'

myped = ped.Ped(ped_file)
myped.addTestFile(field='ind_id', file_pat=test_set_pat)
myped.ped.dropna(subset=['test'], inplace=True)
myped.ped.reset_index(inplace=True)
clm_p = myped.getAllProbands()

# known vars
kv = pandas.read_table(kv_file)
kv['var_id'] = kv.ind_id + '_' +\
               kv.chr + '_' +\
               kv.pos.astype(str)
print kv.groupby(['status', 'descr']).apply(len)
kv.groupby(['status', 'descr']).apply(len).to_csv('/mnt/scratch/asalomatov/data/columbia/kv_sum0.csv')
kv_vcf = pandas.read_table(kv_vcf_file)
kv_vcf = kv_vcf[['ind_id','CHROM', 'POS', 'REF_offspring', 'ALT_base_offspring', 'status', 'descr', 'DP_offspring', 'DP_father', 'DP_mother']]
kv_vcf['var_id'] = kv_vcf.ind_id + '_' +\
               kv_vcf.CHROM + '_' +\
               kv_vcf.POS.astype(str)
kv['in_vcf'] = kv.var_id.isin(kv_vcf.var_id).astype(int)
kv = kv[kv.ind_id.isin(clm_p)]
kv.groupby(['in_vcf', 'descr']).apply(len)
                        


