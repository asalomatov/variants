import os, sys, pandas
from variants import func


batchN = sys.argv[1]

b = pandas.read_csv(
    '/mnt/scratch/asalomatov/data/SPARK/info/ids_by_batch.csv',
    dtype=str)
print b.head()
print b.shape
b_new = pandas.read_table(
    '/mnt/scratch/asalomatov/data/SPARK/info/sampleID_lookup_table_b5',
    dtype=str)
b_new.columns = ['lab_id', 'bn', 'x']
b_new['batch'] = batchN
b_new['bn'] = b_new.bn + '.realigned.recal.bam'
b_new['bam'] = '/mnt/ceph/spark//BCM_trios_BATCH5/' + b_new.bn 
labIS2spID = func.readYml(
    '/mnt/scratch/asalomatov/data/SPARK/info/baylor_id2sp_descr.yml')
b_new['SP_id'] = b_new.lab_id.apply(lambda i: labIS2spID[i])
print b_new.head()
b_new.drop('x', inplace=True, axis=1)
b = pandas.concat([b, b_new])


