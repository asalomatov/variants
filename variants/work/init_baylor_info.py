import os, sys, pandas
#from variants import func

b = pandas.read_csv(
    '/mnt/xfs1/scratch/asalomatov/data/SPARK/info/ids_by_batch_deduplicated.csv',
    dtype=str)
labID2spID = {}
spID2labID = {}
lab2batch = {}
lab2bn = {}
bam2spID = {}

def spId2labId(sp_id, sp2lab_dict):
    if sp_id == '0':
        return sp_id
    elif sp_id[:2] == 'SP':
        return sp2lab_dict[sp_id]
    else:
        return None

for i, row in b.iterrows():
    labID2spID[row['lab_id']] = row['SP_id']
    spID2labID[row['SP_id']] = row['lab_id']
    lab2batch[row['lab_id']] = row['batch']
    lab2bn[row['lab_id']] = row['bn']
    bam2spID[row['bn']] = row['SP_id']
