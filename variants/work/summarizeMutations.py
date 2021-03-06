# from __future__ import print_function
import sys, os
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import pandas, numpy
import func, ped, train, variants
import yaml
import datetime

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
                      u'DP_offspring',
                      u'DP_father',
                      u'DP_mother',
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
                      u'spidex_dpsi_max_tissue',
                      u'spidex_dpsi_zscore',
                      u'spidex_gene',
                      u'spidex_strand',
                      u'spidex_transcript',
                      u'spidex_exon_number',
                      u'spidex_location',
                      u'spidex_cds_type',
                      u'spidex_ss_dist',
                      u'FILTER',
                      u'var_id',
                      u'v_id']


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

def summarizeMutations(infile,
                       outp_dir,
                       config_file,
                       exac_anno='/mnt/scratch/asalomatov/data/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'):
    with open(config_file, 'r') as f:
        cfg = yaml.safe_load(f)

    ped_file = cfg['ped_file']
    myped = ped.Ped(cfg['ped_file'])
    myped.addVcf(file_pat=cfg['vcf_pattern'])
    myped.ped.dropna(subset=['vcf'], inplace=True)
    myped.ped.reset_index(inplace=True)

    #kv_vcf = pandas.read_csv('/mnt/scratch/asalomatov/data/columbia/feature_sets/known/all_known.txt', sep='\t')
    #kv_vcf = kv_vcf[['ind_id','CHROM', 'POS', 'REF_offspring', 'ALT_base_offspring', 'status', 'descr', 'DP_offspring', 'DP_father', 'DP_mother']]
    #kv_vcf = kv_vcf[kv_vcf.descr.isin(['after'])]
    #kv_vcf['var_id'] = kv_vcf.ind_id.astype(str)+'_'+kv_vcf.CHROM.astype(str)+'_'+kv_vcf.POS.astype(str)


    #effects_of_interest = effects_loss_of_func + '|' + effect_damaging_missense + '|' + effect_synon


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
        '|'.join(cfg['snpeff']['effect_synon'])), 'effect_cat'] = 'syn' 
    vn.ix[vn['ANN[*].EFFECT'].str.contains(
        '|'.join(cfg['snpeff']['effect_dmgmis'])), 'effect_cat'] = 'mis' 
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

    vn_full = vn
    allele_frac = vn.alt_DP.astype(float)/vn.DP
    vn = vn[(allele_frac > cfg['alt_allele_frac_range'][0]) & (allele_frac < cfg['alt_allele_frac_range'][1])]

    print '\nallele fraction, pred_labels value_counts:'
    print vn.pred_labels.value_counts()
    print 'allel fraction vars, test_labels value_counts:'
    print vn.status.value_counts()
    calcMetr(vn, msg='all fraction metrics')
    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='allele_frac')])

    vn = vn.replace('ZZZ', '.')
    vn['FILTER'] = vn.apply(func.getFieldFromVCF, args=(myped,), axis=1)
    vn = vn[~vn.FILTER.isnull()]



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
    #vn_FN = vn_FN[~vn_FN.v_id.duplicated()]

    c_TP = (vn.pred_labels == 1) & vn.status.isin(['Y'])
    vn_TP = vn[cols_to_output][c_TP]
    #vn_TP = vn_TP[~vn_TP.v_id.duplicated()]
    
    var_type = cfg['variant_type']
    outp_suffix = '{:%Y-%m-%d_%H-%M-%S-%f}'.format(datetime.datetime.now())

    def writeVariants(df, cols_to_output, var_type, prefix, suffix, outp_dir):
        if df.empty:
            print('%s is empty' % prefix)
            return None
        df[cols_to_output].to_csv(os.path.join(outp_dir,
                                                    '_'.join([prefix,
                                                              var_type,
                                                              suffix,
                                                              '.csv'])), index=False)
    writeVariants(vn, cols_to_output[:-2], var_type, 'ALL', outp_suffix, outp_dir)
    writeVariants(vn_FN, cols_to_output[:-2], var_type, 'FN', outp_suffix, outp_dir)
    writeVariants(vn_TP, cols_to_output[:-2], var_type, 'TP', outp_suffix, outp_dir)
    writeVariants(vn_mis, cols_to_output[:-2], var_type, 'MIS', outp_suffix, outp_dir)
    writeVariants(vn_lof, cols_to_output[:-2], var_type, 'LOF', outp_suffix, outp_dir)
    writeVariants(vn_syn, cols_to_output[:-2], var_type, 'SYN', outp_suffix, outp_dir)
    writeVariants(vn_diff, cols_to_output[:-2]+['step'], var_type, 'DIFF', outp_suffix, outp_dir)
        
#    vn_TP[cols_to_output[:-2]].to_csv(
#        os.path.join(outp_dir, 'true_pos_snp' + outp_suffix + '.csv'), index=False)
#    vn_mis[cols_to_output[:-2]][c_new].to_csv(
#        os.path.join(outp_dir, 'dmg_missense' + outp_suffix + '.csv'), index=False)
#    vn_lof[cols_to_output[:-2]][c_new].to_csv(
#        os.path.join(outp_dir, 'lof' + outp_suffix + '.csv'), index=False)
#    vn_syn[cols_to_output[:-2]][c_new].to_csv(
#        os.path.join(outp_dir, 'syn' + outp_suffix + '.csv'), index=False)
#    vn_diff[cols_to_output[:-2] + ['step']].to_csv(
#        os.path.join(outp_dir, 'lostTP' + outp_suffix + '.csv'), index=False)

    cfg['predictions_file'] = infile

    with open(os.path.join(outp_dir,'cfg' + outp_suffix + '.yml'), 'w') as f:
        yaml.dump(cfg, f, default_flow_style=False)
    

if __name__ == '__main__':
    infile = sys.argv[1]
    outp_dir = sys.argv[2]
    config_file = sys.argv[3]
    func.runInShell('mkdir -p ' + outp_dir)
    summarizeMutations(infile,
                       outp_dir,
                       config_file)





    #cfg = {'population_AF': 0.01,
    #       'snpeff': {'effect_lof': ['exon_loss_variant',
    #                                 'frameshift_variant',
    #                                 'stop_gained',
    #                                 'stop_lost',
    #                                 'start_lost',
    #                                 'splice_acceptor_variant'
    #                                 'splice_donor_variant',
    #                                 'splice_region_variant'],
    #                  'effect_dmgmis': ['missense_variant'],
    #                  'effect_synon': ['synonymous_variant'],
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
    #                },
    #       'min_DP': 10,
    #       'vcf_pattern': '/mnt/scratch/asalomatov/data/columbia/vcf/%s_%s-02_%s-01.annotated-norm.vcf.gz',
    #       'bam_pattern': '/mnt/scratch/asalomatov/data/columbia/bam/Sample.%s.bam',
    #       'ped_file': '/mnt/scratch/asalomatov/data/columbia/pcgc_ped.txt',
    #       'known_variants': '/mnt/scratch/asalomatov/data/columbia/pcgc_denovo_snp.tsv',
    #       'output_directory': '/mnt/scratch/asalomatov/data/columbia/feature_sets_01'
    #}


