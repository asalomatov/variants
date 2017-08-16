import glob, os
import pandas

clms = ['CHROM', 'POS', 'REF', 'ALT', 'effect', 'gene']
x = glob.glob('/mnt/ceph/spark/BCM_trios_BATCH10/*.hq')
df_l = []

for i in x:
    lab_id = os.path.basename(i).split('.')[0]
    print lab_id
    if os.path.getsize(i) > 0:
        df = pandas.read_table(i, header=None, usecols=range(6))
        df.columns = clms
        df['lab_id'] = lab_id
        df_l.append(df)

df = pandas.concat(df_l)

