import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import variants
import ped
import func
import features_vcf as fv
import pandas as pd
import features_vcf
import numpy as np
import os

variants = reload(variants)
func = reload(func)
ped = reload(ped)

###SSC
infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
myped = ped.Ped(infile_ped, ['collection'])
myped.getParents(11006)
#myped.addVcf()
myped.addVcf(file_pat = '/mnt/scratch/asalomatov/data/SSC/vcf/raw/%s.family.vqsr.sorted.vcf.gz')
myped.ped.head()

###Columbia
infile_ped = '/mnt/scratch/asalomatov/data/columbia/pcgc_ped.txt'
myped = ped.Ped(infile_ped, [])
myped.getParents('1-00034')
myped.getFather('1-00034')
myped.
myped.ped.head()
myped.ped.shape
myped.addVcf(file_pat = '/mnt/scratch/asalomatov/data/columbia/vcf/deco/%s_%s-02_%s-01.annotated-deco.vcf.gz')
sum(myped.ped.vcf.notnull())

infile_vcf = '/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/rerun200fam/11006-HC-pm50-ann.vcf.gz'
#infile_vcf = '/mnt/scratch/asalomatov/data/columbia/vcf/deco/1-03173_1-03173-02_1-03173-01.annotated-deco.vcf.gz'
myvars = variants.Variants(infile_vcf, '11006')
myvars.readFromVcf()



myvars.samples
record = myvars.vcf_reader.next()
record.samples
myvars._colNamesFormat()
myvars.describeFormatFields()


###features
###extract all variants from rerun to be tested
infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
myped = ped.Ped(infile_ped, ['collection'])
myped.getParents(11006)
myped.addVcf()
fv = reload(features_vcf)
train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_all_snp.txt"
ftrs = fv.FeaturesVcf(myped, train_set) 
ftrs.variants.head()
df_l = []
for f in fam_trio:
    print f
    fv = reload(features_vcf)
    train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_all_snp.txt"
    ftrs = fv.FeaturesVcf(myped, train_set) 
    if not os.path.isfile('/mnt/scratch/asalomatov/data/SSC/vcf/raw/' + str(f) +'.family.vqsr.sorted.vcf.gz'):
        continue
    df = ftrs._fileHasVariant('/mnt/scratch/asalomatov/data/SSC/vcf/raw/' + str(f) +'.family.vqsr.sorted.vcf.gz', fam_id=f, ind_id=str(f)+'.p1', chrom=None, pos_start=None, pos_end=None)
    df_l.append(df)

len(df_l)
df_test = pd.concat(df_l)
type(df_test)
df_test.isnull().sum()
df_test.shape
df_test.offspring_gt_type.value_counts()
c1 = df_test.vartype.isin(['snp'])
sum(c1)
c3 = df_test['offspring_gt_type'] == 0
df_test = df_test[c1 & ~c3]
df_test.shape
df_test = addVar(df_test, 'ind_id')
ftrs.variants = addVar(ftrs.variants, 'ind_id')
df_test = pd.merge(df_test, ftrs.variants[['var', 'status']], how='left', on='var')
df_test.status.value_counts()
ftrs.variants.status.value_counts()
ftrs.variants.head()
df_test = df_test.dropna(subset=['format_father_0_PL'])
df_test = df_test.dropna(subset=['format_father_1_PL'])
df_test = df_test.dropna(subset=['format_father_2_PL'])
df_test.shape
df_test = df_test.dropna(subset=['format_mother_0_PL'])
df_test = df_test.dropna(subset=['format_mother_1_PL'])
df_test = df_test.dropna(subset=['format_mother_2_PL'])
df_test.shape
df_test = df_test.dropna(subset=['format_offspring_0_PL'])
df_test = df_test.dropna(subset=['format_offspring_1_PL'])
df_test = df_test.dropna(subset=['format_offspring_2_PL'])
df_test.shape
df_test.isnull().sum()

a = pd.DataFrame(ftrs.variants.ind_id.str.split('.',1).tolist(), columns=['fam_id', 'memb'], index=ftrs.variants.index)
ftrs.variants = ftrs.variants.join(a)
ftrs.variants.head()
ftrs.variants.dtypes
ftrs.variants['fam_id'] =  ftrs.variants['fam_id'].astype(int)
sum(ftrs.variants.fam_id.isin(fam_trio))

df_test_num = df2sklearn(df_test)
df_test_num.isnull().sum()
df_test_num.head()
df_test_num.dtypes
df_train_set_num = pd.concat([df_train_num, df_neg_num])
df_train_set_num.head()
df_train_set_num.dtypes
df_train_set_num.shape
df_train_set_num.status01.value_counts()
df_train_set_num.isnull().sum()
df_train_set_num = df_train_set_num.dropna()

ftrs.variants[ftrs.variants..status.value_counts()

### extract some negative examples
train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_verified_snp.txt"
train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_all_snp.txt"
ftrs = fv.FeaturesVcf(myped, train_set) 
ftrs.variants.head()
ftrs.variants.status.value_counts()
ftrs.variants.groupby(['descr', 'status']).apply(len)

df_l = []
for f in df_train['family_id'][:5]:
    print f
    train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_verified_snp.txt"
    ftrs = fv.FeaturesVcf(myped, train_set) 
    df = ftrs._fileHasVariant('/mnt/scratch/asalomatov/data/SSC/vcf/raw/' + str(f) +'.family.vqsr.sorted.vcf.gz',
            fam_id=f, ind_id=str(f)+'.p1', chrom=None, pos_start=None, pos_end=None)
    df_l.append(df)

df = pd.concat(df_l)
type(df)
df.isnull().sum()
df.reset_index(drop=True, inplace=True)
df.isnull().sum()
df.head()
df.shape
df.columns
df.FILTER.value_counts()
df_neg = df
#get some negative examples
#low quality snp
c1 = df.vartype.isin(['snp'])
sum(c1)
c2 = df.FILTER.isin(['PASS'])
sum(c1 & ~c2)
df_lo_qual_snp = df[c1 & ~ c2]
df_lo_qual_snp.isnull().sum()
#child is homRef
df.offspring_gt_type.value_counts()
c1 = df.vartype.isin(['snp'])
sum(c1)
c2 = df.FILTER.isin(['PASS'])
c3 = df['offspring_gt_type'] == 0
df_p1_homref = df[c1 & c2 & c3]
df_p1_homref.isnull().sum()
df_p1_homref.reset_index(drop=True, inplace=True)
N = len(df_p1_homref.index) - 1
N
N_rows = list(set(np.random.randint(low=0, high=N, size=1000)))
len(N_rows)
len(set(N_rows))
df_p1_homref = df_p1_homref.ix[N_rows]
df_p1_homref.shape
df_p1_homref.isnull().sum()
###concat negative examples
df_neg = pd.concat([df_lo_qual_snp, df_p1_homref])
df_neg.shape
df_neg.isnull().sum()
df_neg['status'] = 'N'
df_neg.columns

### annotate training set with features from vcf files
train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_verified_snp.txt"
train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_all_snp.txt"
ftrs = fv.FeaturesVcf(myped, train_set) 
ftrs.ped.ped.head()
df_train = ftrs.extractFeatures()
len(df_train)
df_train = pd.concat(df_train)

df_train.shape
df_train.columns
df_train.status

def addVar(df, field):
    try: df['var'] = df[field].map(str) + '_' + df.POS.map(str) + '_' +  df.REF.map(str)  + '_' +  df.ALT.map(str) 
    except: df['var'] = df[field].map(str) + '_' + df.pos.map(str) + '_' +  df.ref.map(str)  + '_' +  df.alt.map(str) 
    return df

df_train = addVar(df_train, 'ind_id')
df_neg = addVar(df_neg, 'ind_id')
ftrs.variants = addVar(ftrs.variants, 'ind_id')
ftrs.variants.status.value_counts()
df_train = pd.merge(df_train, ftrs.variants[['var', 'status']], how='left', on='var')
df_train.status.value_counts()
df_train['status'][df_train['status'] == 'ND'] = 'Y'
df_neg = pd.merge(df_neg, ftrs.variants[['var', 'status']], how='left', on='var')
df_neg.status.value_counts()
df_neg.shape
df_neg = df_neg[df_neg.status.isnull()]
df_neg.shape
df_neg.status = 'N'
df_train.tail()
df_train.shape

fam_trio = [11193, 11195, 11198, 11827, 13415, 11989, 13733, 11055, 11056, 11545, 11303, 12073,
        12521, 11660, 11388, 11262, 11707, 13008, 12933, 13844, 11184, 11834, 12437, 12430,
        11109, 12532, 11023, 11375, 13314, 13557, 13158, 12300, 11471, 13494, 13857, 12381,
        11205, 13914, 13757, 12015, 13610, 14292, 12157, 13863, 13678, 11120, 13530, 13532,
        11124, 12641, 11083, 11218, 13668, 13742, 11518, 13741, 13333, 12249, 11009, 11510,
        12086, 12674, 11599, 13031, 11096, 11948, 11093, 11947, 11556, 11346, 11224, 13207,
        12444, 11506, 11504, 12036, 11587, 12237, 12335, 12130, 11425, 12238, 14020, 12621,
        13517, 11753, 12185, 11006, 11069, 11141, 12744, 11064, 11148, 11734, 11863, 12225,
        12341, 12346, 12198, 11526, 11523, 13812, 11480, 11928, 12114, 12118, 11246, 12752,
        12296, 12212, 14006, 11498, 11043, 12555, 12667, 13822, 12603, 11396, 11257, 13701,
        11398, 13274, 11653, 11843, 11969]
fam_quad = [13188, 14011, 11964, 13048, 11491, 13793, 11190, 13890, 13835, 12810, 12390, 13169, 12905, 11569, 11629, 11469, 12106, 11773, 13447, 12161, 13116, 11013, 11872, 11172, 11711, 11715, 12011, 14201, 12741, 11390, 11959, 13926, 13335, 11942, 13815, 12373, 12285, 13593, 12703, 11029, 11659, 11472, 11459, 11610, 11788, 13606, 11229, 13346, 11452, 11479, 11722, 13629, 12152, 12153, 12630, 12578, 11696, 12304, 13533, 12358, 12233, 11691]
len(fam_trio)
df_train.count
df_train = df_train[~df_train.family_id.isin(fam_trio)]
df_train.shape
df_train.head()
df_train.isnull().sum()

col_to_keep = [u'QUAL', u'info_BaseQRankSum', u'info_ClippingRankSum', u'info_DP', u'info_FS', u'info_GC', u'info_HRun',
        u'info_MQ', u'info_MQ0', u'info_MQRankSum', u'info_QD',
        u'info_ReadPosRankSum', u'info_SOR', u'format_father_ref_AD', u'format_father_alt_AD', u'format_father_DP',
        u'format_father_GQ', u'format_father_0_PL', u'format_father_1_PL',
        u'format_father_2_PL',u'format_mother_ref_AD', u'format_mother_alt_AD', u'format_mother_DP',
        u'format_mother_GQ', u'format_mother_0_PL', u'format_mother_1_PL', u'format_mother_2_PL',
        u'format_offspring_ref_AD', u'format_offspring_alt_AD', u'format_offspring_DP', u'format_offspring_GQ',
        u'format_offspring_0_PL', u'format_offspring_1_PL', u'format_offspring_2_PL']
col_to_keep_rerun = [u'QUAL', u'info_BaseQRankSum', u'info_ClippingRankSum', u'info_DP', u'info_FS', 
        u'info_MQ', u'info_MQ0', u'info_MQRankSum', u'info_QD',
        u'info_ReadPosRankSum', u'format_father_ref_AD', u'format_father_alt_AD', u'format_father_DP',
        u'format_father_GQ', u'format_father_0_PL', u'format_father_1_PL',
        u'format_father_2_PL',u'format_mother_ref_AD', u'format_mother_alt_AD', u'format_mother_DP',
        u'format_mother_GQ', u'format_mother_0_PL', u'format_mother_1_PL', u'format_mother_2_PL',
        u'format_offspring_ref_AD', u'format_offspring_alt_AD', u'format_offspring_DP', u'format_offspring_GQ',
        u'format_offspring_0_PL', u'format_offspring_1_PL', u'format_offspring_2_PL']

len(col_to_keep_rerun)
len(col_to_keep)
len(feature_cols)

def df2sklearn(mydf, col_to_keep):
    if 'status' in mydf.columns:
        mydf['status01'] = 1
        mydf['status01'][mydf['status'] == 'N'] = 0
        col_to_keep += ['status01']
    col_to_keep = list(set(col_to_keep).intersection(set(mydf.columns)))
    print col_to_keep
    #res = mydf[col_to_keep]
    mydf[col_to_keep] = mydf[col_to_keep].astype(float)
    mydf = mydf.dropna(subset = col_to_keep)
    return mydf[col_tokeep] 

df_train_num = func.df2sklearn(df_train,col_to_keep_rerun)
df_train_num.isnull().sum()
df_train_num.head()
df_train_num.dtypes
df_train_num.status01.value_counts()

df_neg.isnull().sum()
df_neg_num = func.df2sklearn(df_neg, col_to_keep_rerun)
df_neg_num.isnull().sum()
df_neg_num.head()
df_neg_num.dtypes
df_neg_num.status01.value_counts()
df_train_set_num = pd.concat([df_train_num, df_neg_num])
df_train_set_num.head()
df_train_set_num.dtypes
df_train_set_num.shape
df_train_set_num.status01.value_counts()
df_train_set_num.isnull().sum()
df_train_set_num = df_train_set_num.dropna()
df_train_set_num.isnull().sum()
df_train_set_num.shape
df_train_set_num.describe()
#df_train_set_num.to_csv("ssc_snp_training_set_no187.csv", sep="\t", header=True, index=False) 
df_train_set_num = pd.read_table("ssc_snp_training_set_no187.csv", sep="\t") 
df_train_set_num.dtypes
#### below is for my rerun
df_train_set_num = df_train_set_num[col_to_keep_rerun]

def addAlleleBalance(mydf):
    mydf['offspring_allele_balance'] = mydf['format_offspring_alt_AD']/(mydf['format_offspring_alt_AD'] + mydf['format_offspring_ref_AD'])
    mydf['father_allele_balance'] = mydf['format_father_alt_AD']/(mydf['format_father_alt_AD'] + mydf['format_father_ref_AD'])
    mydf['mother_allele_balance'] = mydf['format_mother_alt_AD']/(mydf['format_mother_alt_AD'] + mydf['format_mother_ref_AD'])
    mydf = mydf.dropna(subset = ['offspring_allele_balance', 'father_allele_balance', 'mother_allele_balance'])
    return mydf

df_train_set_num = addAlleleBalance(df_train_set_num)
#df_train_set_num.to_csv("ssc_snp_training_set.csv", sep="\t", header=True, index=False) 

feature_cols = [x for x in df_train_set_num.columns if x not in 'status01']
print feature_cols
response_col = 'status01'
##add allel ballance

### sklearn

#from sklearn.datasets import make_hastie_10_2
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor



# generate synthetic data from ESLII - Example 10.2
#X, y = make_hastie_10_2(n_samples=5000)
#X_train, X_test, y_train, y_test = train_test_split(X, y)
#
# fit estimator
est = GradientBoostingClassifier(n_estimators=2000, max_depth=1, learning_rate=.01)
est_regr = GradientBoostingRegressor(n_estimators=2000, max_depth=1, learning_rate=.01)
len(feature_cols)
X_train = df_train_set_num[feature_cols].values
y_train = df_train_set_num[response_col].values
X_tr, X_te, y_tr, y_te = train_test_split(X_train, y_train)
X_tr
est.fit(X_tr, y_tr)
est_regr.fit(X_tr, y_tr)
feature_importance = pd.DataFrame({'contrib': est.feature_importances_ ,'name': feature_cols})
feature_importance = pd.DataFrame({'contrib': est_regr.feature_importances_ ,'name': feature_cols})
feature_importance.sort(['contrib'], ascending=[False], inplace=True)
feature_importance
#feature_importance.to_excel('feature_contrib.xls', index=False)
est.loss_
sum(est.feature_importances_)

# predict class labels
pred = est.predict(X_te)
pred_regr = est_regr.predict(X_te)
len(pred)
pred_ser = pd.Series(pred)
pred_ser.describe()
pred_ser.value_counts()
sum(y_te == 1)
sum(pred_regr > .17)
sum(pred[pred_regr > .17] == 1)

# score on test data (accuracy)
acc = est.score(X_te, y_te)
print('ACC: %.4f' % acc)

# predict class probabilities
est.predict_proba(X_te)[0]


df = None
df_list = []
len(fam_trio)
#for f in fam_trio:
#for f in fam_quad:
for f in fam_trio + fam_quad:
    memb = 'p1'
    print f
#    if str(f) in files_missing:
#        continue
    fv = reload(features_vcf)
    train_set = "/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_exome_all_snp.txt"
    ftrs = fv.FeaturesVcf(myped, train_set) 
    #filepath ='/mnt/scratch/asalomatov/data/SSC/vcf/raw/' + str(f) +'.family.vqsr.sorted.vcf.gz'
    #filepath ='/mnt/ceph/asalomatov/SSC_Eichler/rerun/ssc' + str(f) + '/' + str(f)  +'-JHC-vars.vcf.gz'
    filepath = '/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/rerun200fam/'+ str(f) + '-JHC-pm50.vcf.gz'
    print filepath
    if not os.path.isfile(filepath):
        continue
    df = ftrs._fileHasVariant(filepath, fam_id=f, ind_id=str(f)+'.'+memb, chrom=None, pos_start=None, pos_end=None)
    if df is None:
        continue
    df_num = func.df2sklearn(df, col_to_keep_rerun)
    df_num = addAlleleBalance(df_num)
    df_num['score'] = est.predict(df_num[feature_cols].values)
#    print 'df_num dim :', df_num.shape
#    df_num = df_num[df_num['score'] > .00]
#    print 'df_num dim :', df_num.shape
    df_list.append(df_num)

len(df_list)
len(df_list_krumm)
type(df)
#df_list_krumm = df_list
df_re = pd.concat(df_list)
df_re.shape
df_re.head()
df_re.score.value_counts()
ftrs.variants.head()
ftrs.variants.shape
c1 = ftrs.variants.ind_id.isin(df_re.ind_id)
sum(c1)
df_kn = ftrs.variants[c1]
df_kn.shape
df_re['var'] = df_re['ind_id'].map(str)+ '_'  + df_re['CHROM'].map(str) + '_' + df_re.POS.map(str) + '_' +  df_re.REF.map(str)  + '_' +  df_re.ALT.map(str) 
ftrs.variants['var'] = ftrs.variants['ind_id'].map(str)+ '_'  + ftrs.variants['chr'].map(str) + '_' + ftrs.variants.pos.map(str) + '_' +  ftrs.variants.ref.map(str)  + '_' +  ftrs.variants.alt.map(str) 

df_re['pos_var'] = df_re['CHROM'].map(str) + '_' + df_re.POS.map(str) + '_' +  df_re.REF.map(str)  + '_' +  df_re.ALT.map(str) 
ftrs.variants['var'] = ftrs.variants['ind_id'].map(str)+ '_'  + ftrs.variants['chr'].map(str) + '_' + ftrs.variants.pos.map(str) + '_' +  ftrs.variants.ref.map(str)  + '_' +  ftrs.variants.alt.map(str) 

df_re.groupby('pos_var').apply(sum)
df_re[df_re.score > .75].groupby('vartype').apply(varsumm)
#how many unique variants?
df_re.vartype.value_counts()
df_re.pheno.value_counts()

len(set(df_re['CHROM'].map(str) + '_' + df_re.POS.map(str) + '_' +  df_re.REF.map(str)  + '_' +  df_re.ALT.map(str)))
c1 = (df_re['dbsnp'] == 'no') & (df_re['vartype'] == 'snp') 
len(set((df_re['CHROM'].map(str) + '_' + df_re.POS.map(str) + '_' +  df_re.REF.map(str)  + '_' +
    df_re.ALT.map(str))[c1]))
###apply hard filter(Ash)
cAsh_snp_1 = df_re['QUAL'] >= 30
cAsh_snp_2 = df_re['info_FS'] < 25
cAsh_snp_3 = df_re['info_QD'] >= 4
cAsh_indel_1 = df_re['info_QD'] >= 1

df_re['status'][cAsh_snp_1 & cAsh_snp_2 & cAsh_snp_3].value_counts()
df_re['dbsnp'][cAsh_snp_1 & cAsh_snp_2 & cAsh_snp_3].value_counts()
df_re['dbsnp'].value_counts()


df_re = pd.merge(df_re, df_kn[['var', 'status', 'descr']], how='left', on='var')
df_re.head()
df_re['status'][df_re.status.isnull()] = 'extra'
df_re['dbsnp'] = 'yes'
df_re['dbsnp'][df_re.ID.isnull()] = 'no'
df_re.ID.isnull().sum()
df_re.dbsnp.value_counts()
df_re.dbsnp.value_counts()
df_re.status.value_counts()
df_re[df_re.score > .75].status.value_counts()
df_re[df_re.score > .75].status.value_counts()
df_re[~df_re.status.isin(['extra'])].dbsnp.value_counts()
df_re.
df_re.groupby('vartype').apply(varsumm)
df_re[df_re.score > .75].groupby('vartype').apply(varsumm)
df_re[df_re.score > .9].groupby('vartype').apply(varsumm)
df_re[df_re.score > .95].groupby('vartype').apply(varsumm)
df_re[df_re.score > .99].groupby('vartype').apply(varsumm)
df_re[df_re.score > .999].groupby('vartype').apply(varsumm)

def varsumm(x):
    validated = sum(x['status'] == 'Y')
    not_determined = sum(x['status'] == 'ND')
    failed = sum(x['status'] == 'N')
    extra = sum(x['dbsnp'] == 'no')
    return pd.Series ([validated, not_determined, failed, extra], index = ['validated', 'not_determined', 'failed', 'extra']) 



df.head()
df.columns
df.shape
df.isnull().sum()
set(df.columns).difference(set(col_to_keep))
set(col_to_keep).difference(set(df.columns))
df_num = df2sklearn(df, col_to_keep_rerun)
df_num.head()
df_num.dtypes
df_num[col_to_keep_rerun].dtypes
df_num = addAlleleBalance(df_num)
df_num = df_num.dropna()
df_num.shape
df.offspring_gt_type.value_counts()
df_num.isnull().sum()
set(df_num.columns).difference(set(df_train_set_num.columns))
set(df_train_set_num.columns).difference(set(df_num.columns))
pred_test = est.predict(df_num[feature_cols].values)
len(pred_test)
df_num['score'] = pred_test
len(pred_test)
pred_test[:100]
pred_test_ser = pd.Series(pred_test)
pred_test_ser.describe()

ftrs.variants[ftrs.variants.ind_id.isin([f])]

