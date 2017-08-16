ssc_dict = {}
wgs40_ssc_dict = {}
with open('/mnt/scratch/asalomatov/data/SSCped/40families.list', 'r') as f:
    for l in f:
        print l
        l_spl = l.strip().split()
        l_spl_ssc = l_spl[0].split('.')
        ssc_dict[l_spl[1]] = (l_spl_ssc[0], l_spl_ssc[1])
        wgs40_ssc_dict[l_spl[0]] = l_spl[1]
    
print ssc_dict
print ssc_dict['SSC05498']
print ssc_dict['SSC05498'][0]
print wgs40_ssc_dict['14590.p1']



import os
import pandas as pd
pd.ExcelFile

pd.algos

pd.E

pd.DataFrame()



eichler_wgs = pd.read_csv( \
        os.path.join(\
        '/mnt/scratch/asalomatov/data/SSCdeNovoCalls/Eichlerlab_Simons_Genome_Project_denovo_calls_06172015.csv'))
eichler_wgs
eichler_wgs.Validation.value_counts()
eichler_wgs.bool
eichler_wgs.abs
eichler_wgs.

eichler_wgs.describe()
eichler_wgs._combine_series_infer
os.
eichler_wgs.rename(columns
eichler_wgs.columns.replace(' ','_', inplace=True)
tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects_of_interest)] 
tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'FAM', 'CHILD', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE']



### Summarize variants from Wendy
pcgc13_full = pd.read_csv('/mnt/scratch/asalomatov/data/columbia/PCGCexomevalidation before March 2013.csv')
pcgc13_full.columns
pcgc13=pcgc13_full[['chromosome', 'position', 'reference_bases', 'alternate_bases', 'blinded_id', 'result', 'denovo']]
pcgc13.head()
pcgc13.shape
pcgc13.columns = ['chr', 'pos', 'ref', 'alt', 'id', 'result' , 'denovo']
pcgc13['source'] = 'pcgc13'

#pcgc13['var'] = pcgc13.chr.map(str) + '_' + pcgc13.pos.map(str) + '_' +  pcgc13.ref.map(str)  + '_' +  pcgc13.alt.map(str) 
#pcgc13['smplvar'] = pcgc13.id.map(str) + '_' + pcgc13.chr.map(str) + '_' + pcgc13.pos.map(str) + '_' +  pcgc13.ref.map(str)  + '_' +  pcgc13.alt.map(str) 
#pcgc13['chr_pos'] = pcgc13.chr.map(str) + '_' + pcgc13.pos.map(str) 
#pcgc13['smpl_chr_pos'] = pcgc13.id.map(str) + '_' + pcgc13.chr.map(str) + '_' + pcgc13.pos.map(str)
#pcgc13_set = set(pcgc13['var'].values)
#pcgc13_id_set = set(pcgc13['id'].values)
#len(pcgc13_id_set)
#pcgc13['result'].value_counts()
#pcgc13['denovo'].value_counts()
#len(pcgc13_id_set)

pcgc14_full = pd.read_csv('/mnt/scratch/asalomatov/data/columbia/PCGC_Variant_Curation_annotation_results_group_summary_11-15-14.csv')
pcgc14_full.columns
pcgc14 = pcgc14_full[['CHR','Coordinates','Reference allele','Alternate allele', 'Blinded ID', 'De novo (yes/no)']]
pcgc14.head()
pcgc14.tail()
pcgc14.dtypes
pcgc14.shape
pcgc14.CHR.value_counts()
pcgc14.columns = ['chr', 'pos', 'ref', 'alt', 'id', 'denovo']
pcgc14['result'] = 'confirmed'
pcgc14['source'] = 'pcgc14'
pcgc14['denovo'].value_counts()

# concat
pcgc = pcgc13
pcgc = pcgc.append(pcgc14)
pcgc.shape
pcgc = pcgc[~(pcgc.ref.isin(['-']) & pcgc.alt.isin(['-']))]
pcgc['pos'] = pcgc['pos'].astype(int)
pcgc.head()
pcgc.dtypes

#convert pos for pcgc13 to hg19 and rerun
#pcgc13_hg19.shape
#pcgc13_hg19.head()
#pcgc13_hg19.columns = ['pos_hg19', 'smplvar0']
#pcgc13_hg19.dtypes
#pcgc['smplvar0'] = pcgc.id.map(str) + '_' + pcgc.chr.map(str) + '_' + pcgc.pos.map(str) + '_' + \
#        pcgc.ref.map(str)  + '_' +  pcgc.alt.map(str) 
#pcgc = pd.merge(pcgc, pcgc13_hg19, how='left', on='smplvar0')
#pcgc.head()
#pcgc.tail()
#pcgc['pos_hg19'][pcgc['pos_hg19'].isnull()] = pcgc[pcgc['pos_hg19'].isnull()]['pos']
#sum(pcgc['pos_hg19'][pcgc['pos_hg19'].isnull()])
#pcgc['pos'] = pcgc['pos_hg19']
#pcgc['pos'] = pcgc['pos'].astype(int)

# extract ref allels and adjust pos for indels
import pysam

def refAtPos(chrom, pos, genref='/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'):
    ref_allel = pysam.faidx(genref, str(chrom)+':'+str(pos)+'-'+str(pos))[1].strip()
    return ref_allel
refAtPos(2, 12340909)

import func
def procIndelRows(row, pos_field):
    if row['ref'] == '-':
        pos = row[pos_field] - 1
        ref_allel = func.refAtPos(row['chr'], pos)
        return pd.Series([pos, ref_allel, ref_allel + row['alt'], 'ins'], index=['pos1', 'ref1', 'alt1', 'vartype'])
    elif row['alt'] == '-':
        pos = row[pos_field] - 1
        ref_allel = func.refAtPos(row['chr'], pos)
        return pd.Series([pos, ref_allel + row['ref'], ref_allel, 'del'], index=['pos1', 'ref1', 'alt1', 'vartype'])
    else:
        pos = row[pos_field]
        return pd.Series([pos, row['ref'], row['alt'], 'snp'], index=['pos1', 'ref1', 'alt1', 'vartype'])

#pcgc.drop(['pos1', 'ref1', 'alt1','vartype', 'smpl_chr_pos'], axis=1, inplace=True)
pcgc.dtypes
x_2 = pcgc.apply(lambda row: procIndelRows(row, 'pos'), axis=1)
type(x_2)
x_2.head()
sum(x_2['pos1'] - pcgc['pos'])
pcgc.head()
pcgc.reset_index(drop=True, inplace=True)
x_2.reset_index(drop=True, inplace=True)
pcgc = pcgc.join(x_2)
pcgc.head()
pcgc.groupby('vartype').apply(len)
#pcgc['chr_pos'] = pcgc.chr.map(str) + '_' + pcgc.pos.map(str) 
#pcgc['var'] = pcgc.chr.map(str) + '_' + pcgc.pos.map(str) + '_' +  pcgc.ref.map(str)  + '_' +  pcgc.alt.map(str) 
#pcgc['smplvar'] = pcgc.id.map(str) + '_' + pcgc.chr.map(str) + '_' + pcgc.pos.map(str) + '_' +  pcgc.ref.map(str)  + '_' +  pcgc.alt.map(str) 
#pcgc['chr_pos'] = pcgc.chr.map(str) + '_' + pcgc.pos.map(str) 
pcgc['smpl_chr_pos'] = pcgc.id.map(str) + '_' + pcgc.chr.map(str) + '_' + pcgc.pos1.map(str)
pcgc['smplvar'] = pcgc.id.map(str) + '_' + pcgc.chr.map(str) + '_' + pcgc.pos1.map(str) + '_' + \
        pcgc.ref1.map(str)  + '_' +  pcgc.alt1.map(str) 
pcgc.shape
pcgc = pcgc[~(pcgc.ref.isin(['-']) & pcgc.alt.isin(['-']))]
pcgc_id_set13 = set(pcgc['id'][pcgc.source.isin(['pcgc13'])].values)
pcgc_id_set14 = set(pcgc['id'][pcgc.source.isin(['pcgc14'])].values)
len(pcgc_id_set13)
len(pcgc_id_set14)
len(pcgc_id_set14.intersection(pcgc_id_set13))
len(pcgc_id_set14.union(pcgc_id_set13))
pcgc[pcgc.vartype.isin(['ins'])]
pcgc[pcgc.vartype.isin(['del'])]

# remove duplicates
pcgc = pcgc.drop_duplicates(subset='smplvar', take_last=True)
pcgc.shape
pcgc['has_file'] = pcgc['id'].isin(vcf_smpl)
pcgc['has_file'].value_counts()
#pcgc['file'] = pcgc['id'] + '_' + pcgc['id'] + '-02_' + pcgc['id'] + '-01.annotated-deco.vcf.gz'
pcgc['file'] = pcgc['id'] + '_' + pcgc['id'] + '-02_' + pcgc['id'] + '-01.annotated.vcf.gz'
pcgc['file'][~pcgc['has_file']] = None
sum(pcgc['file'].isnull())

pcgc['file_variant'] = None
sum(pcgc['file'].isnull())
type(pcgc.id)
pcgc.id.head()
pcgc.head()

### now extract variants from corresponding vcf files

#for i, row in pcgc.iterrows():
#    if not row['file'] is None:
#        fname = '/mnt/scratch/asalomatov/data/columbia/vcf/deco/' + row['file']
#        chrom = row['chr']
#        pos_end = row['pos']
#        pos_start = pos_end - 1
#        print chrom, pos_start, pos_end, fname
#        myvars = Variants(fname, row['id'], chrom, pos_start, pos_end)
#        if not next(myvars.vcf_reader, None) is None:
#            row['file_has_variant'] = True


x_1 = pcgc[pcgc['has_file']==True].apply(lambda row: func.fileHasVariant('/mnt/scratch/asalomatov/data/columbia/vcf/' + row['file'], row['id'], \
        row['chr'], row['pos1'] - 1, row['pos1']), axis=1)

type(x_1)
x_1.head()
x_1.tail()
len(x_1)
print(x_1)
[x for x in x_1.tolist() if len(x)>1]
[x for x in x_1.tolist() if '00' in x]
import itertools
len( list(itertools.chain.from_iterable(x_1.tolist())))
x_1_merged = list(itertools.chain.from_iterable(x_1))
x_1_merged[:10]
[x for x in x_1_merged if '00' in x]
len(x_1_merged)
x_1_merged

pcgc['file_has_variant'] = pcgc['smplvar'].isin(x_1_merged)
pcgc['file_has_variant'].value_counts()
pcgc_0['file_has_variant'].value_counts()
pcgc.head()
pcgc.tail()
pcgc.dtypes
c1 = (~pcgc.file_has_variant) & (pcgc['file_has_variant'] != pcgc['chr_pos'])
c2 = (pcgc.file_has_variant.notnull()) & (pcgc['file_has_variant'] == pcgc['chr_pos'])
pcgc['variant_found'] = False
pcgc['variant_found'][c2] = True
pcgc['variant_found'].value_counts()
sum(c1)
sum(c2)
pcgc[(pcgc.file_has_variant.notnull()) & (pcgc['file_has_variant'] != pcgc['var'])].tail()
pcgc[(pcgc['file_has_variant'].notnull())]
type(pcgc.iloc[612]['file_has_variant'])
sum(pcgc.file_has_variant == pcgc.var)

def mySum(x):
    total = len(x['id'])
    return pd.Series ([total], index = ['total'])

pcgc_sum = pcgc[pcgc.has_file].groupby(['source', 'vartype', 'result', 'file_has_variant']).apply(mySum)
pcgc_sum.to_excel
pcgc_sum.to_excel('pcgc_sum.xls', index=False)
pcgc.to_excel('pcgc.xls', index=False)
pcgc_0.groupby(['source', 'vartype', 'result', 'has_file', 'file_has_variant']).apply(mySum)
pcgc_0.groupby(['source', 'result', 'has_file', 'file_has_variant']).apply(mySum)
pcgc_0[pcgc_0.has_file].groupby(['source', 'result', 'has_file', 'file_has_variant']).apply(mySum)
pcgc_0[pcgc_0.has_file].groupby(['source', 'vartype', 'result', 'file_has_variant']).apply(mySum)
sum(pcgc.file_has_variant.isnull())
pcgc_0.shape
pcgc.shape
pcgc.head()
pcgc_0.head()



#### 

pcgc[pcgc.has_file & pcgc.vartype.isin(['ins', 'del']) & ~pcgc.file_has_variant]
pcgc_0 = pcgc

# create bed file for pcgc13 to liftOver from hg18 import to hg19
y = pcgc[pcgc['source']=='pcgc13'][['chr','pos', 'smplvar']]
y
y['pos_st'] = y['pos'] -1 
y[['chr', 'pos_st', 'pos']]
y['chr'] = 'chr' + y['chr'] 

y[['chr', 'pos_st', 'pos', 'smplvar']].to_csv("pcgc13.bed", sep="\t", header=False, index=False)

### read pcgc13 converted to hg19 from hg18 import with liftOver
pcgc13_hg19 = pd.read_table('hglft_genome_4393_d68f00.bed', header=None)
pcgc13_hg19.head()
pcgc13_hg19.columns = ['chr', 'pos', 'pos_hg19', 'smplvar']
pcgc13_hg19 = pcgc13_hg19[[ 'pos_hg19', 'smplvar']]

# now convert pcgc to hg19 by substituting pcgc13 coordinates
pcgc['pos_hg19'] = pcgc['pos1']
pcgc.head()
pcgc.drop('pos_hg19', axis=1, inplace=True)
y = pd.merge(pcgc, pcgc13_hg19, how='left', on='smplvar')
y.head()
y.shape
y[y.source == 'pcgc14']['pos_hg19'].head()
y[y.source == 'pcgc14']['pos_hg19'].tail()
y['pos_hg19'][y.source == 'pcgc14'] = y[y.source == 'pcgc14']['pos']
y['pos_hg19'][y['pos_hg19'].isnull()] = y[y['pos_hg19'].isnull()]['pos']
y['pos_hg19'] = y['pos_hg19'].astype(int) 

x = y[y['has_file']==True].apply(lambda row: fileHasVariant('/mnt/scratch/asalomatov/data/columbia/vcf/deco/' + row['file'], row['id'], \
        row['chr'], row['pos_hg19'] - 1, row['pos_hg19']), axis=1)

x
y.head()
y.tail()
y.dtypes
y['file_has_variant'] = None
y['file_has_variant'][y['has_file']==True] = x
y['chr_pos_hg19'] = y['chr'].astype(str) + '_' + y['pos_hg19'].astype(str) 
c1 = (y.file_has_variant.notnull()) & (y['file_has_variant'] != y['chr_pos_hg19'])
c2 = (y.file_has_variant.notnull()) & (y['file_has_variant'] == y['chr_pos_hg19'])
y['variant_found'] = False
y['variant_found'][c2] = True
y['variant_found'].value_counts()

y.head()

def mySum(x):
    total = len(x['file_has_variant'])
    return pd.Series ([total], index = ['total'])

mysum = pcgc.groupby(['source', 'result', 'has_file', 'variant_found']).apply(mySum)
mysum_hg19 = y.groupby(['source', 'result', 'has_file', 'variant_found']).apply(mySum)


y[(y['source']=='pcgc13') & (y['variant_found']) & (y['result']=='not_confirmed')]
pcgc[(pcgc['source']=='pcgc13') & (pcgc['variant_found']) & (pcgc['result']=='not_confirmed')]
y[(y['id']=='1-00174')]

y[(y['source']=='pcgc13') & (y['file_has_variant'].notnull()) & (y['result']=='not_confirmed')]
pcgc[(pcgc['source']=='pcgc13') & (pcgc['file_has_variant']True) & (pcgc['result']=='not_confirmed')]
y.dtypes



### process Iioss and Krumm variants
ssc_dnv = pd.read_table('/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/SSC_exome_denovo_Ios_Krumm.csv',
        sep="\t", index_col=False, skiprows=1, header=None)
ssc_dnv.columns
ssc_dnv.shape
ssc_dnv.dtypes
ssc_dnv.tail()
ssc_dnv.Study.value_counts()
ssc_dnv = ssc_dnv[[4,0,1,2,3,11,12]]
ssc_dnv.columns = ['ind_id', 'chr', 'pos', 'ref', 'alt', 'status', 'descr']
ssc_dnv.head()
ssc_dnv['vartype'] = 'snp'
x = ssc_dnv['alt'].apply(len)
y = ssc_dnv['ref'].apply(len)
ssc_dnv['vartype'][(x > 1)]= 'del'
ssc_dnv['vartype'][(y > 1)]= 'ins'
ssc_dnv['vartype'].value_counts()
ssc_dnv['status'].value_counts()
ssc_dnv[ssc_dnv['vartype'] == \
    'snp'].to_csv("/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_all_snp.txt", sep="\t", header=True, index=False)

ssc_dnv_known_status = ssc_dnv[ssc_dnv['status'] != 'ND']
ssc_dnv_known_status.shape
ssc_dnv_known_status.shape
ssc_dnv_known_status.head()
ssc_dnv_known_status.tail()
ssc_dnv_known_status['vartype'].value_counts()
ssc_dnv_known_status[ssc_dnv_known_status['vartype'] != 'snp'].to
ssc_dnv_known_status[ssc_dnv_known_status['vartype'] !=
        'snp'].to_csv("/mnt/scratch/asalomatov/data/SSCdeNovoCalls/ssc_exome_verified_indels.txt", sep="\t", header=True, index=False)
ssc_dnv_known_status[ssc_dnv_known_status['vartype'] == \
    'snp'].to_csv("/mnt/scratch/asalomatov/data/SSCdeNovoCalls/ssc_exome_verified_snp.txt", sep="\t", header=True, index=False)
ssc_dnv_known_status.groupby(['vartype', 'status']).apply(len)



### dbSNP fields
myfile = '/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz'
myvars = variants.Variants(myfile, 'aaa')
myvars.describeInfoFields()


