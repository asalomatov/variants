import pandas
import os
import sys


def vipPed(row, i):
    if row['sex'] == 'male':
        sex = '1'
    else:
        sex = '2'
    if row['relationship_to_iip'] == 'iip':
        pheno = '2'
    else:
        pheno = '1'

    res = pandas.DataFrame({
        'fam_id':
        [
            str(i),
            str(i),
            str(i)],
        'ind_id':
        [
            row['individual'],
            row['father'],
            row['mother']],
        'fa_id':
        [
            row['father'],
            '0',
            '0'],
        'mo_id':
        [
            row['mother'],
            '0',
            '0'],
        'sex':
        [
            sex,
            '1',
            '2'],
        'pheno':
        [
            pheno,
            '1',
            '1']})
    print res
    return res

sys.exit('done')

b=pandas.read_csv('/mnt/xfs1/home/asalomatov/baylor_sample_output_file_batch'+str(btch)+'.csv')
b['ind_id'] = b.Internal_Procesing_SampleID
# b['fa_id'] = 0
# b['mo_id'] = 0

b['batch_num'] = btch


def genTrioName(n):
    if n < 10:
        return 'trio0000' + str(n)
    elif n >= 10 and n < 100:
        return 'trio000' + str(n)
    elif n >= 100 and n < 1000:
        return 'trio00' + str(n)
    elif n >= 1000 and n < 10000:
        return 'trio0' + str(n)
    else:
        return 'trio' + str(n)


def getBamName(row):
    bam_name = os.path.join(bam_dir,
                            'batch_' + str(row['batch_num']),
                            row['File_name']+'.realigned.recal.bam')
    bai_name = bam_name + '.bai'
    res = pandas.Series([bam_name, bai_name], ['bam', 'bai'], dtype=str)
    return res


b_d = {}
for i, row in b.iterrows():
    b_d[row['SPARK_Child_ID']] = row['ind_id']

fa_id = []
for i, row in b.iterrows():
    if row['SPARK_BioFather_ID'] in b_d:
        fa_id.append(b_d[row['SPARK_BioFather_ID']])
    else:
        fa_id.append(0)
b['fa_id'] = fa_id

mo_id = []
for i, row in b.iterrows():
    if row['SPARK_BioMother_ID'] in b_d:
        mo_id.append(b_d[row['SPARK_BioMother_ID']])
    else:
        mo_id.append(0)
b['mo_id'] = mo_id

b['sex'] = 1
b.ix[b.SPARK_Child_Sex == 'F', 'sex'] = 2

b['pheno'] = 1
b.ix[b.SPARK_Child_Phenotype, 'pheno'] = 2

trio_dict = {}
for i, row in b.iterrows():
    if row['mo_id'] != 0 and row['fa_id'] != 0:
        tr_n = prev_tr_num + 1
        trio_dict[row['ind_id']] = genTrioName(tr_n)
        trio_dict[row['fa_id']] = genTrioName(tr_n)
        trio_dict[row['mo_id']] = genTrioName(tr_n)
        prev_tr_num = tr_n

trio_num = []
for i, row in b.iterrows():
    trio_num.append(trio_dict[row['ind_id']])

b['fam_id'] = trio_num

b = b.merge(b.apply(getBamName, axis=1), left_index=True, right_index=True)

b['vcf_hc'] = None
b.ix[b.batch_num == 1, 'vcf_hc'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch1/hc/C8-HC-vars-flr.vcf.gz'
b.ix[b.batch_num == 2, 'vcf_hc'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch2/hc/C9-HC-vars-flr.vcf.gz'
b['vcf_jhc'] = None
b.ix[b.batch_num == 1, 'vcf_jhc'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch1/jhc/C8-JHC-vars.vcf.gz'
b.ix[b.batch_num == 2, 'vcf_jhc'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch1/jhc/C9-JHC-vars.vcf.gz'
b['vcf_fb'] = None
b.ix[b.batch_num == 1, 'vcf_fb'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch1/fb/C8-FB-vars.vcf.gz'
b.ix[b.batch_num == 2, 'vcf_fb'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch1/fb/C9-FB-vars.vcf.gz'
b['vcf_pl'] = None
b.ix[b.batch_num == 1, 'vcf_pl'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch1/pl/C8-PL-vars.vcf.gz'
b.ix[b.batch_num == 2, 'vcf_pl'] = '/mnt/xfs1/home/asalomatov/projects/spark/batch1/pl/C9-PL-vars.vcf.gz'

b = b[['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno', 'bam', 'bai', 'vcf_hc', 'vcf_jhc', 'vcf_fb', 'vcf_pl']]
b.sort_values(by='fam_id', inplace=True)
cols =['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno', 'bam']
b[cols + ['vcf_hc']].to_csv(
    '/mnt/scratch/asalomatov/data/SPARK/spark_hc_b' + str(btch) + '.ped',
    sep='\t', header=False, index=False)
b[cols + ['vcf_jhc']].to_csv(
    '/mnt/scratch/asalomatov/data/SPARK/spark_jhc_b' + str(btch) + '.ped',
    sep='\t', header=False, index=False)
b[cols + ['vcf_pl']].to_csv(
    '/mnt/scratch/asalomatov/data/SPARK/spark_pl_b' + str(btch) + '.ped',
    sep='\t', header=False, index=False)
b[cols + ['vcf_fb']].to_csv(
    '/mnt/scratch/asalomatov/data/SPARK/spark_fb_b' + str(btch) + '.ped',
    sep='\t', header=False, index=False)


'''

x.columns=['ind', 'bam_id', 'smpl', 'x', 'y']
x
    
x['bam']=x.apply(lambda x: runInShell('ls /mnt/ceph/spark/BCM_trios/'+x+'*.bam', True))
x['bam']=x.bam_id.apply(lambda x: runInShell('ls /mnt/ceph/spark/BCM_trios/'+x+'*.bam', True))
x
x['sf_id']=None
x
x.ix[x.smpl.isin([641943]), 'bam']
x.ix[x.smpl.isin([641943]), 'sf_id'] = 'trio001.p1_641943'
x
x.ix[x.smpl.isin([641964]), 'sf_id'] = 'trio001.fa_641964'
x.ix[x.smpl.isin([641942]), 'sf_id'] = 'trio001.mo_641942'
x
x.ix[x.smpl.isin([641945]), 'sf_id'] = 'trio002.p1_641945'
x.ix[x.smpl.isin([641946]), 'sf_id'] = 'trio002.fa_641946'
x.ix[x.smpl.isin([641947]), 'sf_id'] = 'trio002.mo_641947'
x
x.ix[x.smpl.isin([642940]), 'sf_id'] = 'trio003.p1_642940'
x.ix[x.smpl.isin([642732]), 'sf_id'] = 'trio003.fa_642732'
x.ix[x.smpl.isin([642941]), 'sf_id'] = 'trio003.mo_642941'
x
x.columns
x.columns[-1]
x.columns[:-1]
x.columns = x.columns[:-1]+['ind_id']
len(x.columns)
len(x.columns[:-1])
len(x.columns[:-1]+['ind_id'])
len(list(x.columns[:-1])+['ind_id'])
(list(x.columns[:-1])+['ind_id'])
x.columns=(list(x.columns[:-1])+['ind_id'])
x
x['fam_id']=x.ind_id.apply(lambda x: x.split('.')[0])
x
get_ipython().magic(u'pwd ')
for i, row in x.iterrows():
    print i , row
    
for i, row in x.iterrows():
    print row.bam
    
for i, row in x.iterrows():
    print row.bam, row.ind_i
    d
    
for i, row in x.iterrows():
    print row.bam, row.ind_id
    
for i, row in x.iterrows():
    runInShell('cd /mnt/ceph/spark/BCM_trios/; ln -s '+row.bam+' '+row.ind_id+'.bam')
    
for i, row in x.iterrows():
    runInShell('cd ~projects/spark/bam; ln -s '+row.bam+' '+row.ind_id+'.bam')
    
for i, row in x.iterrows():
    runInShell('cd ~/projects/spark/bam; ln -s '+row.bam+' '+row.ind_id+'.bam')
    
for i, row in x.iterrows():
    runInShell('cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s '+row.bam+' '+row.ind_id+'.bam')
    
for i, row in x.iterrows():
    print 'cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s '+row.bam+' '+row.ind_id+'.bam'
    
for i, row in x.iterrows():
    runInShell('cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s '+row.bam+' '+row.ind_id+'.bam')
    
get_ipython().magic(u'll /mnt/scratch/asalomatov/data/SPARK/bam/')
get_ipython().system(u'rm  /mnt/scratch/asalomatov/data/SPARK/bam/*.bam')
for i, row in x.iterrows():
    runInShell(' 'join(['cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s', row.bam, row.ind_id+'.bam'])
    
    
    )
    
for i, row in x.iterrows():
    runInShell(' '.join(['cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s', row.bam, row.ind_id+'.bam']))
    
get_ipython().system(u'rm  /mnt/scratch/asalomatov/data/SPARK/bam/*.bam')
for i, row in x.iterrows():
    print (' '.join(['cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s', row.bam, row.ind_id+'.bam']))
    
get_ipython().system(u'cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s /mnt/ceph/spark/BCM_trios/C8CU5ACXX-2-ID01.realigned.recal.bam')
trio001.mo_641942.bam
for i, row in x.iterrows():
    print (' '.join(['cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s', row.bam, row.ind_id+'.bam']))
    
for i, row in x.iterrows():
    print '\n' in  (' '.join(['cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s', row.bam, row.ind_id+'.bam']))
    
x.bam.apply(lambda x: '\n' in x)
x.bam.str.replace('\n', '')
x.ix [:, 'bam'] = x.bam.str.replace('\n', '')
x.ind_id.apply(lambda x: '\n' in x)
for i, row in x.iterrows():
    runInShell(' '.join(['cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s', row.bam, row.ind_id+'.bam']))
    
get_ipython().magic(u'll /mnt/scratch/asalomatov/data/SPARK/bam/')
get_ipython().magic(u'rm /mnt/scratch/asalomatov/data/SPARK/bam/*.bam')
get_ipython().magic(u'll /mnt/scratch/asalomatov/data/SPARK/bam/')
for i, row in x.iterrows():
    runInShell(' '.join(['cd /mnt/scratch/asalomatov/data/SPARK/bam/; ln -s', row.bam, row.ind_id+'.bam']))
    
get_ipython().magic(u'll /mnt/scratch/asalomatov/data/SPARK/bam/')
x
x[x.fam_id=='trio001']
x[x.fam_id=='trio002']
x[x.fam_id=='trio003']
get_ipython().system(u'head /mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped')
x.head
x.head()
x['fa_id']='0'
x.ind_id.str.contains('p1')
isProband = x.ind_id.str.contains('p1')
x['pheno']=1
x.ix[isProband, 'pheno'] = 2
def asFfaForFam(fam_id, df);
def asFfaForFam(fam_id, df):
    x = df.ind_id[df.fam_id==fam_id]
    return x
asFfaForFam('trio001', x)
def asFfaForFam(fam_id, df):
    x = df.ind_id[df.fam_id==fam_id]
    return x[x.str.contains('fa')]
asFfaForFam('trio001', x)
type(asFfaForFam('trio001', x))
def asFfaForFam(fam_id, df):
    x = df.ind_id[df.fam_id==fam_id]
    return str(x[x.str.contains('fa')])
type(asFfaForFam('trio001', x))
asFfaForFam('trio001', x)
def asFfaForFam(fam_id, df):
    x = df.ind_id[df.fam_id==fam_id]
    return list(x[x.str.contains('fa')])[0]
asFfaForFam('trio001', x)
def asFmoForFam(fam_id, df):
    x = df.ind_id[df.fam_id==fam_id]
    return list(x[x.str.contains('mo')])[0]
asFmoForFam('trio001', x)
asFmoForFam('trio002', x)
x.columns
get_ipython().magic(u'ls ')
x['mo_id']='0'
x[x.ind_id.str.contains('p1'),'mo_id']=x.fam_id[x.ind_id.str.contains('p1')0.apply(kaka)

)
def asFmoForFam(fam_id, df=x):
    x = df.ind_id[df.fam_id==fam_id]
    return list(x[x.str.contains('mo')])[0]
def asFfaForFam(fam_id, df=x):
    x = df.ind_id[df.fam_id==fam_id]
    return list(x[x.str.contains('fa')])[0]
x.fam_id[x.ind_id.str.contains('p1')].apply(asFfaForFam)
x[x.ind_id.str.contains('p1'),'fa_id']
x.ix(x.ind_id.str.contains('p1'),'fa_id')
x.ix[x.ind_id.str.contains('p1'),'fa_id']
x.ix[x.ind_id.str.contains('p1'),'fa_id'] = x.fam_id[x.ind_id.str.contains('p1')].apply(asFfaForFam)
x.ix[x.ind_id.str.contains('p1'),'mo_id'] = x.fam_id[x.ind_id.str.contains('p1')].apply(asFmoForFam)
x
x['sex']=1
x
x['sex']=2
x.ix[x.ind_id.str.contains('p1') | x.ind_id.str.contains(, 'sex']
)
x.ix[x.ind_id.str.contains('p1') | x.ind_id.str.contains('fa'), 'sex'] = 1
x
x[['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno']]
x.sort_values('ind_id')
x
x = x.sort_values('ind_id')
x[['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno']]
get_ipython().set_next_input(u"x[['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno']].to_csv");get_ipython().magic(u'pinfo to_csv')
x[['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno']].to_csv
get_ipython().magic(u'pinfo pandas.DataFrame.to_csv')
x[['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno']].to_csv('/mnt/scratch/asalomatov/data/SPARK/spark.ped', header=False, index=False, sep='\t')
import yaml
cfg=yaml.load('cfg_spark.yml')
cfg
f=open('cfg_spark.yml', 'r')
cfg=yaml.load(f)
cfg=yaml.load(f)
cfg
f
cfg=yaml.load(f)
print cfg
f.close()
f=open('cfg_spark.yml', 'r')
cfg=yaml.load(f)
cfg
cfg['known_variants'] is None

'''
