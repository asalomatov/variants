import pandas
import numpy
import yaml
import os
import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants/work/')
from init_baylor_info import spID2labID
from init_baylor_info import labID2spID
from init_baylor_info import lab2batch
from init_baylor_info import lab2bn


def trioIDforInd(df, ind_id):
    res = df['Trio_Descriptor'][df.SP_Descriptor == ind_id]
    return res.iloc[0]


def phenoForInd(df, ind_id, trio_d, trio_d_clm='spark_trio_id'):
    res = df['person_affected'][(df.SP_Descriptor == ind_id) &
                                (df[trio_d_clm].str.contains(trio_d))]
    res = res.iloc[0]
    if str(res) == 'True':
        return '2'
    else:
        return '1'


def sexForInd(df, ind_id, trio_d, trio_d_clm='spark_trio_id'):
    res = df['person_sex'][(df.SP_Descriptor == ind_id) &
                           (df[trio_d_clm].str.contains(trio_d))]
    res = res.iloc[0]
    if res == 'M':
        return '1'
    elif res == 'F':
        return '2'
    else:
        return '0'


def roleIDForTrio(df, trio_d, role_d, trio_d_clm='spark_trio_id'):
    print trio_d
    res = df['SP_Descriptor'][(df[trio_d_clm].str.contains(trio_d)) &
                              (df.role.apply(lambda i: i.lower()) == role_d)]
    print res
    assert(len(set(res)) == 1)
    return res.iloc[0]


def retKthOrFirst(y, k):
    y_spl = str(y).split(',')
    if len(y_spl) == 1:
        return y
    else:
        return y_spl[k]


def pedRow(row, x):
    print 'pedRow'
    print row['Mnemonic_Descriptor']
    print row['role']
    mnem = row['Mnemonic_Descriptor'].split(',')
    for k, m in enumerate(mnem):
        spark_id = row['spark_trio_id'].split(',')[k]
        if row['role'].lower() == 'child':
            res = pandas.Series([m[:6],
                                 row['SP_Descriptor'],
                                 roleIDForTrio(x, spark_id, 'father'),
                                 roleIDForTrio(x, spark_id, 'mother'),
                                 sexForInd(x, row['SP_Descriptor'],
                                           spark_id),
                                 phenoForInd(x, row['SP_Descriptor'],
                                             spark_id),
                                 row['lab_sample_id']],
                                ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                                 'sex', 'pheno', 'lab_sample_id'])
        if row['role'].lower() in ['father', 'mother']:
            res = pandas.Series([m[:6],
                                 row['SP_Descriptor'],
                                 '0',
                                 '0',
                                 sexForInd(x, row['SP_Descriptor'],
                                           spark_id),
                                 phenoForInd(x, row['SP_Descriptor'],
                                             spark_id),
                                 row['lab_sample_id']],
                                ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                                 'sex', 'pheno', 'lab_sample_id'])
    return res


def sexForSPID(spid, df):
    res = list(set(df.person_sex[df.SP_Descriptor == spid]))
    assert(len(res) == 1)
    if res[0] == 'F':
        return '2'
    elif res[0] == 'M':
        return '1'
    else:
        return '0'


def labIDForSPID(spid, df):
    return None
#    res = list(set(df.lab_sample_id[df.SP_Descriptor == spid]))
#    assert(len(res) == 1)
#    return(res[0])


def phenoForSPID(spid, df):
    res = list(set(df.person_affected[df.SP_Descriptor == spid]))
    assert(len(res) == 1)
    if str(res[0]) == 'False':
        return '1'
    elif str(res[0]) == 'True':
        return '2'
    else:
        return '0'


def pedFaForChild(row, x):
    res = pandas.Series([row['fam_id'],
                         row['fa_id'],
                         '0',
                         '0',
                         sexForSPID(row['fa_id'], x),
                         phenoForSPID(row['fa_id'], x),
                         labIDForSPID(row['fa_id'], x)],
                        ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                         'sex', 'pheno', 'lab_sample_id'])
    return res


def pedMoForChild(row, x):
    res = pandas.Series([row['fam_id'],
                         row['mo_id'],
                         '0',
                         '0',
                         sexForSPID(row['mo_id'], x),
                         phenoForSPID(row['mo_id'], x),
                         labIDForSPID(row['mo_id'], x)],
                        ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                         'sex', 'pheno', 'lab_sample_id'])
    return res


def pedRow_trioID(row):
    if row['role'].lower() == 'child':
        res = pandas.Series([row['spark_trio_id'],
                             row['SP_Descriptor'],
                             roleIDForTrio(x, row['spark_trio_id'],
                                           'father', 'spark_trio_id'),
                             roleIDForTrio(x, row['spark_trio_id'],
                                           'mother', 'spark_trio_id'),
                             sexForInd(x, row['SP_Descriptor'],
                                       row['spark_trio_id'], 'spark_trio_id'),
                             phenoForInd(x, row['SP_Descriptor'],
                                         row['spark_trio_id'], 'spark_trio_id'),
                             row['lab_sample_id']],
                            ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                             'sex', 'pheno', 'lab_sample_id'])
    if row['role'].lower() in ['father', 'mother']:
        res = pandas.Series([row['spark_trio_id'],
                             row['SP_Descriptor'],
                             '0',
                             '0',
                             sexForInd(x, row['SP_Descriptor'],
                                       row['spark_trio_id'], 'spark_trio_id'),
                             phenoForInd(x, row['SP_Descriptor'],
                                         row['spark_trio_id'], 'spark_trio_id'),
                             row['lab_sample_id']],
                            ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                             'sex', 'pheno', 'lab_sample_id'])
    return res


def insertBamAndBatch(row, id_column='lab_sample_id', N_batches=3):
    bam_path = '/mnt/scratch/asalomatov/data/SPARK/bam'
    with open('/mnt/scratch/asalomatov/data/SPARK/info/baylor_id2bam.yml', 'r') as f:
        bay2bam = yaml.safe_load(f)
    b_num = numpy.arange(N_batches) + 1
    b_num = b_num.astype(str)
    try:
        bam_name = bay2bam[int(row[id_column])] + '.realigned.recal.bam'
        for b in b_num:
            abs_bam = os.path.join(bam_path, 'batch_' + b, bam_name)
            if os.path.isfile(abs_bam):
                return pandas.Series(['batch_' + b, abs_bam], ['batch', 'bam'])
    except Exception:
        sys.exc_clear()
    return pandas.Series([None, None], ['batch', 'bam'])




### change below

x = pandas.read_csv(
    '/mnt/xfs1/scratch/asalomatov/data/SPARK/info/trio_samples_all_menemonic_20170620.csv',
    dtype=str)

b = pandas.read_csv(
    '/mnt/xfs1/scratch/asalomatov/data/SPARK/info/trios_500.csv', header=None)
b.columns = ['SP_id']
#b = pandas.read_csv(
#    '/mnt/xfs1/scratch/asalomatov/data/SPARK/info/ids_by_batch_deduplicated.csv')
#b = pandas.read_table(
#    '/mnt/xfs1/scratch/asalomatov/data/SPARK/info/sampleID_lookup_table_b10', sep='\s')

#b.columns = ['bam', 'smpl_id', 'barcode', 'lab_id']
#x_b = x[x.lab_sample_id.isin(b.lab_id)]
x_b = x[x.SP_Descriptor.isin(b.SP_id)]
x_b = x_b[~x_b.duplicated()]
x_b['fam_id'] = x_b.Mnemonic_Descriptor.apply(lambda i: i[:6])


def roleFromMnem(mnem):
    md = {'c': 'Child', 'm': 'Mother', 'f': 'Father'}
    return md[mnem[-1:]]

x_b['role'] = x_b.Mnemonic_Descriptor.apply(roleFromMnem)

ch = x_b[x_b.role == 'Child']
pch = ch.apply(lambda row: pedRow(row, x_b), axis=1)
pfa = pch.apply(lambda row: pedFaForChild(row, x_b), axis=1)
pmo = pch.apply(lambda row: pedMoForChild(row, x_b), axis=1)
p = pandas.concat([pch, pfa, pmo])
p.sort_values(['fam_id', 'ind_id'], inplace=True)
p['bam'] = None  # p.lab_sample_id.apply(lambda i: lab2bn[i])
p.columns = [u'fam_id', u'ind_id', u'fa_id', u'mo_id', u'sex', u'pheno',
       u'lab_id', u'bam']
sys.exit('stop!')




