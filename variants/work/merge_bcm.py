import os
import sys
import pandas
import glob
from variants import  func
from init_baylor_info import spID2labID
from init_baylor_info import labID2spID
from init_baylor_info import lab2batch

inp_dir = sys.argv[1:]
d_l = []
for i_d in inp_dir:
    print 'dir is ', i_d
    f_l = glob.glob(os.path.join(i_d, '*.denovo.hq'))
    print 'file list is ', f_l
    for f in f_l:
        try:
            d = pandas.read_table(f, header=None, dtype=str, usecols=range(6))
        except:
            continue
        d.columns = ['CHROM', 'POS', 'REF', 'ALT', 'effect', 'gene']
        d['lab_id'] = os.path.basename(f).split('.')[0]
        d_l.append(d)

df = pandas.concat(d_l,)
df.reset_index(inplace=True, drop=True)
#labID2spID = func.readYml('/mnt/scratch/asalomatov/data/SPARK/info/baylor_id2sp_descr.yml')
df['SP_id'] = df.lab_id.apply(lambda z: labID2spID[z])
df[['lab_id', 'SP_id'] + ['CHROM', 'POS', 'REF', 'ALT', 'effect', 'gene']].to_csv(
    '/mnt/ceph/users/asalomatov/spark/bcm/bcm_denovo_hq_b1-6.csv',
    index=False)
