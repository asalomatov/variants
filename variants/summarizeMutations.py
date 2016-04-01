#!python
import pandas, numpy
import sys, os
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import train
import yaml

def calcMetr(vn_df, msg=' '):
    tst = train.TrainTest('x',
                      '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt',
                      ['status'],
                      ['descr'])
    vn_d = vn_df[~vn_df.var_id.duplicated()]
    tst.pred_y = numpy.array(vn_d['pred_labels'].astype(int))
    tst.test_set_y = tst.pred_y * 0
    tst.test_set_y[numpy.array(vn_d.status.isin(['Y']))] = 1
    print '\n\n HHHHHHHAAAAAAA'
    print pandas.Series(tst.test_set_y).value_counts()
    print msg
    tst.getMetrics()

def getDiff(df_full, df_new, msg, field='var_id'):
    df_full_d = df_full[~df_full[field].duplicated()]
    df_new_d = df_new[~df_new[field].duplicated()]
    res = df_full_d[df_full_d.status.isin(['Y']) &
                   df_full_d.pred_labels.isin([1]) &
                   (~df_full_d[field].isin(df_new_d[field]))]
    if not res.empty: res.ix[:, 'step'] = msg
    return res

infile = sys.argv[1]
outp_dir = sys.argv[2]
outp_suffix = sys.argv[3]
exac_anno =\
'/mnt/scratch/asalomatov/data/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'
#cfg = {'population_AF': 0.001,
#       'snpeff': {'effect_lof': effects_loss_of_func.split('|'),
#                  'effect_dmgmis': effect_damaging_missense.split('|'),
#                  'effect_synon': effect_synon.split('|'),
#                  'biotype': ['protein_coding'],
#                  'impact_lof': ['HIGH'],
#                  'impact_dmgmis': [],
#                  'genes': []},
#       'max_cohort_freq': 1,
#       'db_nsfp': { 'metaSVM_pred': ['D'],
#                    'cadd_phred': 25,
#                    'combined': {'cadd_phred': 15,
#                                 'polyphen2_pred': ['D'],
#                                 'sift_pred': ['D']}
#                }
#}
with open('cfg.yml', 'r') as f:
    cfg = yaml.safe_load(f)

with open(os.path.join(outp_dir,'cfg' + outp_suffix + '.yml'), 'w') as f:
    yaml.dump(cfg, f, default_flow_style=False)

#kv_vcf = pandas.read_csv('/mnt/scratch/asalomatov/data/columbia/feature_sets/known/all_known.txt', sep='\t')
#kv_vcf = kv_vcf[['ind_id','CHROM', 'POS', 'REF_offspring', 'ALT_base_offspring', 'status', 'descr', 'DP_offspring', 'DP_father', 'DP_mother']]
#kv_vcf = kv_vcf[kv_vcf.descr.isin(['after'])]
#kv_vcf['var_id'] = kv_vcf.ind_id.astype(str)+'_'+kv_vcf.CHROM.astype(str)+'_'+kv_vcf.POS.astype(str)


#effects_of_interest = effects_loss_of_func + '|' + effect_damaging_missense + '|' + effect_synon
cols_to_output = [u'CHROM',
                  u'POS',
                  u'ID',
                  u'REF',
                  u'ALT',
                  u'ind_id',
                  u'effect_cat',
                  u'pred_labels',
                  u'pred_prob',
                  u'status',
                  u'ref_DP',
                  u'alt_DP',
                  u'DP',
                  u'ANN[*].EFFECT',
                  u'ANN[*].IMPACT',
                  u'ANN[*].GENE',
                  u'ANN[*].GENEID',
                  u'ANN[*].BIOTYPE',
                  u'dbNSFP_rs_dbSNP142',
                  u'dbNSFP_1000Gp3_AF',
                  u'dbNSFP_ExAC_AF',
                  u'dbNSFP_Polyphen2_HVAR_pred',
                  u'dbNSFP_Polyphen2_HDIV_pred',
                  u'dbNSFP_CADD_phred',
                  u'dbNSFP_MetaSVM_pred',
                  u'dbNSFP_SIFT_pred',
                  u'syn_z',
                  u'mis_z',
                  u'lof_z',
                  u'pLI',
                  u'pRec',
                  u'pNull',
                  u'var_id',
                  u'v_id']


exac = pandas.read_table(exac_anno)
vn = pandas.read_table(infile)
vn.columns = vn.columns.str.translate(None, '#')
print vn.shape
vn.ix[:, 'gene'] = vn['ANN[*].GENE']
vn = vn.merge(
    exac[[u'syn_z', u'mis_z', u'lof_z', u'pLI', u'pRec', u'pNull', u'gene']],
    on='gene', how='left')
print vn.shape
#sys.exit(1)
vn['v_id'] = vn.ind_id + '_' +\
             vn['CHROM'].astype(str) + '_' +\
             vn.POS.astype(str) + '_' +\
             vn['ANN[*].GENE'] + '_' +\
             vn['ANN[*].EFFECT'] + '_' +\
             vn['ANN[*].IMPACT']
vn['var_id'] = vn.ind_id + '_' +\
             vn['CHROM'].astype(str) + '_' +\
             vn.POS.astype(str)
vn['chr_pos'] = vn['CHROM'].astype(str) + '_' +\
             vn.POS.astype(str)


#vn = vn.merge(kv_vcf[['var_id', 'status']], on='var_id', how='left')
#print vn.shape

vn = vn[~vn.v_id.duplicated()]
# stats before any filtering
print '\ndeduped and annotated vars, pred_labels value_counts:'
print vn.pred_labels.value_counts()
print 'deduped and annotated vars, test_labels value_counts:'
print vn.status.value_counts()
calcMetr(vn, msg='deduped metrics')
vn_full = vn


var_freq = vn.groupby('chr_pos').apply(lambda x: len(x['ind_id'].unique()))
var_freq_2 = var_freq[var_freq > cfg['max_cohort_freq']]
vn = vn[~vn.chr_pos.isin(var_freq_2.index)]
print '\ncohort freq vars, pred_labels value_counts:'
print vn.pred_labels.value_counts()
print 'cohort freq vars, test_labels value_counts:'
print vn.status.value_counts()
calcMetr(vn, msg='cohort_freq')
vn_diff =  getDiff(vn_full, vn, msg='cohort_freq')



vn.ix[:, 'effect_cat'] = None
vn.ix[vn['ANN[*].EFFECT'].str.contains(
    '|'.join(cfg['snpeff']['effect_dmgmis'])), 'effect_cat'] = 'mis' 
vn.ix[vn['ANN[*].EFFECT'].str.contains(
    '|'.join(cfg['snpeff']['effect_synon'])), 'effect_cat'] = 'syn' 
vn.ix[vn['ANN[*].EFFECT'].str.contains(
    '|'.join(cfg['snpeff']['effect_lof'])), 'effect_cat'] = 'lof' 
print vn.shape
vn_full = vn
vn = vn.dropna(subset=['effect_cat'], axis=0)
print vn.shape
#vn = vn[vn['ANN[*].EFFECT'].str.contains(effects_of_interest)] 

print '\neffects of interest vars, pred_labels value_counts:'
print vn.pred_labels.value_counts()
print 'effects of interest vars, test_labels value_counts:'
print vn.status.value_counts()
calcMetr(vn, msg='effects metrics')
vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='effects')])


vn.ix[vn.dbNSFP_1000Gp3_AF.isin(['.']), 'dbNSFP_1000Gp3_AF'] = '0'
vn.ix[vn.dbNSFP_ExAC_AF.isin(['.']), 'dbNSFP_ExAC_AF'] = '0'
vn.ix[:, 'dbNSFP_1000Gp3_AF'] = vn.dbNSFP_1000Gp3_AF.str.replace('ZZZ', '0')
vn.ix[:, 'dbNSFP_ExAC_AF'] = vn.dbNSFP_ExAC_AF.str.replace('ZZZ', '0')
vn.ix[:, 'dbNSFP_1000Gp3_AF'] = vn.dbNSFP_1000Gp3_AF.str.replace(',.', ',0')
vn.ix[:, 'dbNSFP_ExAC_AF'] = vn.dbNSFP_ExAC_AF.str.replace(',.', ',0')
vn.ix[:, 'dbNSFP_1000Gp3_AF'] = vn.dbNSFP_1000Gp3_AF.str.replace('.,', '0,')
vn.ix[:, 'dbNSFP_ExAC_AF'] = vn.dbNSFP_ExAC_AF.str.replace('.,', '0,')
vn.ix[:, 'dbNSFP_1000Gp3_AF'] = vn.dbNSFP_1000Gp3_AF.apply(lambda x: min(map(float, x.split(','))))
vn.ix[:, 'dbNSFP_ExAC_AF'] = vn.dbNSFP_ExAC_AF.apply(lambda x: min(map(float, x.split(','))))
vn.ix[:, 'dbNSFP_1000Gp3_AF'] = vn.dbNSFP_1000Gp3_AF.astype(float)
vn.ix[:, 'dbNSFP_ExAC_AF'] = vn.dbNSFP_ExAC_AF.astype(float)

vn_full = vn

vn = vn[(vn.dbNSFP_1000Gp3_AF < cfg['population_AF']) &
        (vn.dbNSFP_ExAC_AF < cfg['population_AF'])]
print '\nAF vars, pred_labels value_counts:'
print vn.pred_labels.value_counts()
print 'AF vars, test_labels value_counts:'
print vn.status.value_counts()
calcMetr(vn, msg='AF metrics')
vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='pop_freq')])

vn_full = vn
vn = vn[vn['ANN[*].BIOTYPE'].str.contains('|'.join(cfg['snpeff']['biotype']))]

print '\nprotein coding vars, pred_labels value_counts:'
print vn.pred_labels.value_counts()
print 'protein coding vars, test_labels value_counts:'
print vn.status.value_counts()
calcMetr(vn, msg='protein coding metrics')
vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='protein')])

vn = vn.replace('ZZZ', '.') 

c_missense = vn['effect_cat'] == 'mis'
c_lof = vn['effect_cat'] == 'lof'
c_syn = vn['effect_cat'] == 'syn'

c_metaSVM_D = vn.dbNSFP_MetaSVM_pred.str.contains('|'.join(cfg['db_nsfp']['metaSVM_pred']))
c_metaSVM_null = vn.dbNSFP_MetaSVM_pred.isin(['ZZZ', '.'])

c_cadd_null = vn.dbNSFP_CADD_phred.isin(['ZZZ', '.'])
c_cadd_D = vn.dbNSFP_CADD_phred[~c_cadd_null].apply(
    lambda x: min(map(float, x.split(',')))) >= cfg['db_nsfp']['cadd_phred']
c_cadd_15 = vn.dbNSFP_CADD_phred[~c_cadd_null].apply(
    lambda x: min(map(float, x.split(',')))) >= cfg['db_nsfp']['combined']['cadd_phred']

c_poly_HVAR_null = vn.dbNSFP_Polyphen2_HVAR_pred.isin(['ZZZ', '.'])
c_poly_HDIV_null = vn.dbNSFP_Polyphen2_HVAR_pred.isin(['ZZZ', '.'])
c_poly_HVAR_D = vn.dbNSFP_Polyphen2_HVAR_pred.str.contains(
    '|'.join(cfg['db_nsfp']['combined']['polyphen2_pred']))
c_poly_HDIV_D = vn.dbNSFP_Polyphen2_HVAR_pred.str.contains(
    '|'.join(cfg['db_nsfp']['combined']['polyphen2_pred']))

c_sift_null = vn.dbNSFP_SIFT_pred.isin(['ZZZ', '.'])
c_sift_D = vn.dbNSFP_SIFT_pred.str.contains(
    '|'.join(cfg['db_nsfp']['combined']['sift_pred']))
c_new = (vn.pred_labels == 1) & (~vn.status.isin(['Y']))

c_dmg_miss = c_metaSVM_D | c_cadd_D | ((c_poly_HDIV_D | c_poly_HVAR_D) & c_sift_D & c_cadd_15)
vn_full = vn[c_missense]
vn_mis = vn[c_dmg_miss & c_missense]
vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn_mis, msg='dmg_miss')])

c_impact_lof = vn['ANN[*].IMPACT'].str.contains(
    '|'.join(cfg['snpeff']['impact_lof']))

vn_full = vn[c_lof]
vn_lof = vn[c_lof & c_impact_lof]
vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn_lof, msg='impact_lof')])

vn_syn = vn[c_syn]
print vn.shape

c_FN = (vn.pred_labels == 0) & vn.status.isin(['Y'])
vn_FN = vn[cols_to_output][c_FN]
vn_FN = vn_FN[~vn_FN.var_id.duplicated()]

c_TP = (vn.pred_labels == 1) & vn.status.isin(['Y'])
vn_TP = vn[cols_to_output][c_TP]
vn_TP = vn_TP[~vn_TP.var_id.duplicated()]


vn[cols_to_output[:-2]].to_csv(
    os.path.join(outp_dir, 'all_snp' + outp_suffix + '.csv'), index=False)
vn_FN[cols_to_output[:-2]].to_csv(
    os.path.join(outp_dir, 'false_neg_snp' + outp_suffix + '.csv'), index=False)
vn_TP[cols_to_output[:-2]].to_csv(
    os.path.join(outp_dir, 'true_pos_snp' + outp_suffix + '.csv'), index=False)
vn_mis[cols_to_output[:-2]][c_new].to_csv(
    os.path.join(outp_dir, 'dmg_missense' + outp_suffix + '.csv'), index=False)
vn_lof[cols_to_output[:-2]][c_new].to_csv(
    os.path.join(outp_dir, 'lof' + outp_suffix + '.csv'), index=False)
vn_syn[cols_to_output[:-2]][c_new].to_csv(
    os.path.join(outp_dir, 'syn' + outp_suffix + '.csv'), index=False)
vn_diff[cols_to_output[:-2] + ['step']].to_csv(
    os.path.join(outp_dir, 'lostTP' + outp_suffix + '.csv'), index=False)


sys.exit(1)









fam_quad = ['13188', '14011', '11964', '13048', '11491', '13793', '11190', '13890', '13835', '12810', '12390', '13169', '12905', '11569', '11629', '11469', '12106', '11773', '13447', '12161', '13116', '11013', '11872', '11172', '11711', '11715', '12011', '14201', '12741', '11390', '11959', '13926', '13335', '11942', '13815', '12373', '12285', '13593', '12703', '11029', '11659', '11472', '11459', '11610', '11788', '13606', '11229', '13346', '11452', '11479', '11722', '13629', '12152', '12153', '12630', '12578', '11696', '12304', '13533', '12358', '12233', '11691']

fam_trio = ['11193', '11195', '11198', '11827', '13415', '11989', '13733', '11055', '11056', '11545', '11303', '12073', '12521', '11660', '11388', '11262', '11707', '13008', '12933', '13844', '11184', '11834', '12437', '12430', '11109', '12532', '11023', '11375', '13314', '13557', '13158', '12300', '11471', '13494', '13857', '12381', '11205', '13914', '13757', '12015', '13610', '14292', '12157', '13863', '13678', '11120', '13530', '13532', '11124', '12641', '11083', '11218', '13668', '13742', '11518', '13741', '13333', '12249', '11009', '11510', '12086', '12674', '11599', '13031', '11096', '11948', '11093', '11947', '11556', '11346', '11224', '13207', '12444', '11506', '11504', '12036', '11587', '12237', '12335', '12130', '11425', '12238', '14020', '12621', '13517', '11753', '12185', '11006', '11069', '11141', '12744', '11064', '11148', '11734', '11863', '12225', '12341', '12346', '12198', '11526', '11523', '13812', '11480', '11928', '12114', '12118', '11246', '12752', '12296', '12212', '14006', '11498', '11043', '12555', '12667', '13822', '12603', '11396', '11257', '13701', '11398', '13274', '11653', '11843', '11969']

families = fam_trio + fam_quad
# read data into dataframes for these families

hc_p1 = pd.DataFrame()
hc_s1 = pd.DataFrame()
jhc_p1 = pd.DataFrame()
jhc_s1 = pd.DataFrame()
fb_p1 = pd.DataFrame()
fb_s1 = pd.DataFrame()
pl_p1 = pd.DataFrame()
pl_s1 = pd.DataFrame()
ios = pd.DataFrame()
wgs = pd.DataFrame()
denovo = pd.DataFrame()
genes = 'ANK2|ASH1L|CHD8|GRIN2B|SCN2A|DSCAM|ADNP|DYRK1A|SHANK3|CHD2' 
effects = 'exon_loss_variant|frameshift_variant|stop_gained|stop_lost|start_lost|splice_acceptor_variant|splice_donor_variant|rare_amino_acid_variant|missense_variant|inframe_insertion|disruptive_inframe_insertion|inframe_deletion|disruptive_inframe_deletion|5_prime_UTR_truncation+exon_loss_variant|3_prime_UTR_truncation+exon_loss|splice_region_variant'
effects_loss_of_func = 'exon_loss_variant|frameshift_variant|stop_gained|stop_lost|start_lost|splice_acceptor_variant|splice_donor_variant|splice_region_variant'
effect_damaging_missense='missense_variant'
effects_of_interest = effects_loss_of_func + '|' + effect_damaging_missense
effect1 = 'frameshift_variant'
effect2 = 'start_lost'
effect3 = 'splice_region_variant'
effect4 = 'stop_gained'

p1_dir = '/mnt/ceph/asalomatov/data/SSCexome/rerun200fam/p1'
s1_dir = '/mnt/ceph/asalomatov/data/SSCexome/rerun200fam/s1'
ios_dir = '/mnt/ceph/asalomatov/data/SSCexome/200fam'
denovo_file = '/mnt/ceph/asalomatov/data/SSCexome/rerun200fam/ios_denovo/ios_denovo-pm50-ann.txt'

if os.path.isfile(os.path.join(denovo_file)):
    tempdf = pd.read_table(os.path.join(denovo_file))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
    tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects_of_interest)] 
    tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
    tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'FAM', 'CHILD', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE']
    denovo = denovo.append(tempdf)
    denovo = denovo[denovo.FAM.isin(families)]

for f in families:
    print 'processing familiy ', f
    if os.path.isfile(os.path.join(ios_dir,str(f)+"-ios-pm50-ann.txt")):
        tempdf = pd.read_table(os.path.join(ios_dir,str(f)+"-ios-pm50-ann.txt"))
#        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
#        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'FAM', 'EFFECT', 'IMPACT', 'GENE',
                'FEATURE', 'MISC_INFO']
        ios = ios.append(tempdf)

for f in fam_trio:
    print 'processing trio ', f
    if os.path.isfile(os.path.join(p1_dir,f+"-HC-pm50-ann-dnmfp1.txt")):
        tempdf = pd.read_table(os.path.join(p1_dir,f+"-HC-pm50-ann-dnmfp1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1']
        tempdf['FAM'] = int(f)
        tempdf['DNMFilt_s1'] = 'NA'
        tempdf['GT_s1'] = '.'
        tempdf['FAM_TYPE'] = 'TRIO'
        hc_p1 = hc_p1.append(tempdf)

        tempdf = pd.read_table(os.path.join(p1_dir,f+"-JHC-pm50-ann-dnmfp1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1']
        tempdf['FAM'] = int(f)
        tempdf['DNMFilt_s1'] = 'NA'
        tempdf['GT_s1'] = '.'
        tempdf['FAM_TYPE'] = 'TRIO'
        jhc_p1 = jhc_p1.append(tempdf)

        tempdf = pd.read_table(os.path.join(p1_dir,f+"-FB-pm50-ann-dnmfp1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1']
        tempdf['FAM'] = int(f)
        tempdf['DNMFilt_s1'] = 'NA'
        tempdf['GT_s1'] = '.'
        tempdf['FAM_TYPE'] = 'TRIO'
        fb_p1 = fb_p1.append(tempdf)

        tempdf = pd.read_table(os.path.join(p1_dir,f+"-PL-pm50-ann-dnmfp1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1']
        tempdf['FAM'] = int(f)
        tempdf['DNMFilt_s1'] = 'NA'
        tempdf['GT_s1'] = '.'
        tempdf['FAM_TYPE'] = 'TRIO'
        pl_p1 = pl_p1.append(tempdf)

for f in fam_quad:
    print 'processing quad', f
    if os.path.isfile(os.path.join(s1_dir,f+"-HC-pm50-ann-dnmfp1-dnmfs1.txt")):
        tempdf = pd.read_table(os.path.join(s1_dir,f+"-HC-pm50-ann-dnmfp1-dnmfs1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'DNMFilt_s1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1', 'GT_s1']
        tempdf['FAM'] = int(f)
        tempdf['FAM_TYPE'] = 'QUAD'
        hc_s1 = hc_s1.append(tempdf)

        tempdf = pd.read_table(os.path.join(s1_dir,f+"-JHC-pm50-ann-dnmfp1-dnmfs1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'DNMFilt_s1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1', 'GT_s1']
        tempdf['FAM'] = int(f)
        tempdf['FAM_TYPE'] = 'QUAD'
        jhc_s1 = jhc_s1.append(tempdf)

        tempdf = pd.read_table(os.path.join(s1_dir,f+"-FB-pm50-ann-dnmfp1-dnmfs1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'DNMFilt_s1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1', 'GT_s1']
        tempdf['FAM'] = int(f)
        tempdf['FAM_TYPE'] = 'QUAD'
        fb_s1 = fb_s1.append(tempdf)

        tempdf = pd.read_table(os.path.join(s1_dir,f+"-PL-pm50-ann-dnmfp1-dnmfs1.txt"))
#    c1 = tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    c2 = tempdf['ANN[*].EFFECT'].str.contains(effect2)
#    c3 = tempdf['ANN[*].EFFECT'].str.contains(effect3)
#    c4 = tempdf['ANN[*].EFFECT'].str.contains(effect4) & ~tempdf['ANN[*].EFFECT'].str.contains(effect1)
#    tempdf = tempdf[c1 | c2 | c3 | c4]
        tempdf = tempdf[tempdf['ANN[*].EFFECT'].str.contains(effects)] 
        tempdf = tempdf[tempdf['ANN[*].GENE'].str.contains(genes)] 
        tempdf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'VARTYPE', 'DNMFilt_p1', 'DNMFilt_s1', 'EFFECT', 'IMPACT', 'GENE', 'FEATURE', 'GT_fa', 'GT_mo', 'GT_p1', 'GT_s1']
        tempdf['FAM'] = int(f)
        tempdf['FAM_TYPE'] = 'QUAD'
        pl_s1 = pl_s1.append(tempdf)

hc = pd.DataFrame()
jhc = pd.DataFrame()
fb = pd.DataFrame()
pl = pd.DataFrame()
hc = hc.append(hc_p1)
hc = hc.append(hc_s1)
jhc = jhc.append(jhc_p1)
jhc = jhc.append(jhc_s1)
fb = fb.append(fb_p1)
fb = fb.append(fb_s1)
pl = pl.append(pl_p1)
pl = pl.append(pl_s1)

#mydfs = [hc, jhc, fb, pl]
#for df in mydfs:
hc.DNMFilt_p1 = hc.DNMFilt_p1.replace(['NA'], [None])
hc.DNMFilt_s1 = hc.DNMFilt_s1.replace(['NA'], [None])
hc.DNMFilt_p1 = hc.DNMFilt_p1.astype(float).fillna(0.0)
hc.DNMFilt_s1 = hc.DNMFilt_s1.astype(float).fillna(0.0)
jhc.DNMFilt_p1 = jhc.DNMFilt_p1.replace(['NA'], [None])
jhc.DNMFilt_s1 = jhc.DNMFilt_s1.replace(['NA'], [None])
jhc.DNMFilt_p1 = jhc.DNMFilt_p1.astype(float).fillna(0.0)
jhc.DNMFilt_s1 = jhc.DNMFilt_s1.astype(float).fillna(0.0)
fb.DNMFilt_p1 = fb.DNMFilt_p1.replace(['NA'], [None])
fb.DNMFilt_s1 = fb.DNMFilt_s1.replace(['NA'], [None])
fb.DNMFilt_p1 = fb.DNMFilt_p1.astype(float).fillna(0.0)
fb.DNMFilt_s1 = fb.DNMFilt_s1.astype(float).fillna(0.0)
pl.DNMFilt_p1 = pl.DNMFilt_p1.replace(['NA'], [None])
pl.DNMFilt_s1 = pl.DNMFilt_s1.replace(['NA'], [None])
pl.DNMFilt_p1 = pl.DNMFilt_p1.astype(float).fillna(0.0)
pl.DNMFilt_s1 = pl.DNMFilt_s1.astype(float).fillna(0.0)

#for df in mydfs:
hc = hc[hc['EFFECT'].str.contains(effects_of_interest)] 
jhc = jhc[jhc['EFFECT'].str.contains(effects_of_interest)] 
fb = fb[fb['EFFECT'].str.contains(effects_of_interest)] 
pl = pl[pl['EFFECT'].str.contains(effects_of_interest)] 
denovo = denovo[denovo['EFFECT'].str.contains(effects_of_interest)] 
ios = ios[ios['EFFECT'].str.contains(effects_of_interest)] 



def oneGene(g):
    g_l = g.split(',')
    for gg in g_l:
        if gg in genes:
            return gg
    return None
def mutCategory(mut):
    mut_l = mut.split(',')
    for mm in mut_l:
        if mm in effects_loss_of_func:
            return 'loss_of_function'
        elif mm in effect_damaging_missense:
            return 'missense'
        else:
            return None

#for df in mydfs:
hc['GENE'] = hc['GENE'].apply(oneGene)
hc = hc[~hc['GENE'].isnull()]
jhc['GENE'] = jhc['GENE'].apply(oneGene)
jhc = jhc[~jhc['GENE'].isnull()]
fb['GENE'] = fb['GENE'].apply(oneGene)
fb = fb[~fb['GENE'].isnull()]
pl['GENE'] = pl['GENE'].apply(oneGene)
pl = pl[~pl['GENE'].isnull()]
denovo['GENE'] = denovo['GENE'].apply(oneGene)
denovo = denovo[~denovo['GENE'].isnull()]
ios['GENE'] = ios['GENE'].apply(oneGene)
ios = ios[~ios['GENE'].isnull()]

hc['EFFECT'] = hc['EFFECT'].apply(mutCategory)
hc = hc[~hc['EFFECT'].isnull()]
jhc['EFFECT'] = jhc['EFFECT'].apply(mutCategory)
jhc = jhc[~jhc['EFFECT'].isnull()]
fb['EFFECT'] = fb['EFFECT'].apply(mutCategory)
fb = fb[~fb['EFFECT'].isnull()]
pl['EFFECT'] = pl['EFFECT'].apply(mutCategory)
pl = pl[~pl['EFFECT'].isnull()]
denovo['EFFECT'] = denovo['EFFECT'].apply(mutCategory)
denovo = denovo[~denovo['EFFECT'].isnull()]
ios['EFFECT'] = ios['EFFECT'].apply(mutCategory)
ios = ios[~ios['EFFECT'].isnull()]

denovo['var'] = denovo.CHROM.map(str) + '_' + denovo.POS.map(str) + '_' + denovo.REF + '_' + denovo.ALT
ios['var'] = ios.CHROM.map(str) + '_' + ios.POS.map(str) + '_' + ios.REF + '_' + ios.ALT
hc['var'] = hc.CHROM.map(str) + '_' + hc.POS.map(str) + '_' + hc.REF + '_' + hc.ALT
jhc['var'] = jhc.CHROM.map(str) + '_' + jhc.POS.map(str) + '_' + jhc.REF + '_' + jhc.ALT
fb['var'] = fb.CHROM.map(str) + '_' + fb.POS.map(str) + '_' + fb.REF + '_' + fb.ALT
pl['var'] = pl.CHROM.map(str) + '_' + pl.POS.map(str) + '_' + pl.REF + '_' + pl.ALT

variants = pd.DataFrame()
variants = variants.append(jhc)
variants = variants.append(hc[~hc['var'].isin(variants['var'])])
variants = variants.append(fb[~fb['var'].isin(variants['var'])])
variants = variants.append(pl[~pl['var'].isin(variants['var'])])
variants['var_in_p1'] = (variants['GT_p1'] == '0/1') | (variants['GT_p1'] == '1/1')
variants['var_in_s1'] = (variants['GT_s1'] == '0/1') | (variants['GT_s1'] == '1/1')

def varsumm(x,dnmCutoff, ios_df, denovo_df):
    missense = sum(x['EFFECT'] == 'missense')
    loss_of_func = sum(x['EFFECT'] == 'loss_of_function')
    N_p1 = sum((x['GT_p1'] == '0/1') | (x['GT_p1'] == '1/1'))
    N_s1 = sum((x['GT_s1'] == '0/1') | (x['GT_s1'] == '1/1'))
    denovo_p1 = sum(x['DNMFilt_p1'] > dnmCutoff)
    denovo_s1 = sum(x['DNMFilt_s1'] > dnmCutoff)
    in_ios_denovo = sum(x['var'].isin(denovo_df['var']))
#    in_ios = sum(x['var'].isin(ios_df['var']))
#    not_in_ios = sum(~x['var'].isin(ios_df['var']))
#    ios_not_sf = sum(~ios_df['var'].isin(x['var']))
#    return pd.Series ([loss_of_func, missense, N_p1, N_s1, denovo_p1, denovo_s1, in_ios_denovo, in_ios, not_in_ios, ios_not_sf], index = ['loss_of_func', 'missense', 'N_p1', 'N_s1', 'denovo_p1', 'denovo_s1', 'in_ios_denovo', 'in_ios', 'not_in_ios', 'ios_not_sf'])
    return pd.Series ([loss_of_func, missense, N_p1, N_s1, denovo_p1, denovo_s1, in_ios_denovo], index = ['loss_of_func', 'missense', 'N_p1', 'N_s1', 'denovo_p1', 'denovo_s1', 'in_ios_denovo']) 

def varsumm1(x,dnmCutoff, ios_df, denovo_df):
    total = len(x['EFFECT'])
    p1 = len(x['EFFECT'][x['var_in_p1']])
    s1 = len(x['EFFECT'][x['var_in_s1']])
    parents = len(x['EFFECT'][x['var_in_parents']])
    p1_and_parents = len(x['EFFECT'][(x['var_in_p1']) & x['var_in_parents']])
    s1_and_parents = len(x['EFFECT'][(x['var_in_s1']) & x['var_in_parents']])
    denovo_p1 = sum(x['DNMFilt_p1'][x['var_in_p1']] > dnmCutoff)
    denovo_s1 = sum(x['DNMFilt_s1'][x['var_in_s1']] > dnmCutoff)
    p1_in_ios_denovo = sum(x['var'][x['var_in_p1']].isin(denovo_df['var']))
    s1_in_ios_denovo = sum(x['var'][x['var_in_s1']].isin(denovo_df['var']))
#    in_ios = sum(x['var'].isin(ios_df['var']))
#    not_in_ios = sum(~x['var'].isin(ios_df['var']))
#    ios_not_sf = sum(~ios_df['var'].isin(x['var']))
#    return pd.Series ([loss_of_func, missense, N_p1, N_s1, denovo_p1, denovo_s1, in_ios_denovo, in_ios, not_in_ios, ios_not_sf], index = ['loss_of_func', 'missense', 'N_p1', 'N_s1', 'denovo_p1', 'denovo_s1', 'in_ios_denovo', 'in_ios', 'not_in_ios', 'ios_not_sf'])
    return pd.Series ([total,p1, s1, parents, p1_and_parents, s1_and_parents, denovo_p1, denovo_s1, p1_in_ios_denovo, s1_in_ios_denovo], index = ['total','p1','s1', 'parents', 'p1_and_parents', 's1_and_parents', 'denovo_p1', 'denovo_s1', 'p1_in_ios_denovo', 's1_in_ios_denovo'])

sum_gene = variants.groupby('GENE').apply(varsumm)
sum_gene_effect = variants.groupby(['GENE', 'EFFECT']).apply(varsumm)




