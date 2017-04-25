from __future__ import print_function
import sys
import os
import pandas
import numpy
import func
import ped
import train
import yaml
# import datetime

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
                  u'allele_frac',
                  u'DP_offspring',
                  u'DP_father',
                  u'DP_mother',
                  u'ANN[*].EFFECT',
                  u'ANN[*].IMPACT',
                  u'ANN[*].GENE',
                  u'ANN[*].GENEID',
                  u'ANN[*].FEATURE',
                  u'ANN[*].FEATUREID',
                  u'ANN[*].BIOTYPE',
                  u'Ann[*].RANK',
                  u'ANN[*].HGVS_C',
                  u'ANN[*].HGVS_P',
                  u'ANN[*].CDNA_POS',
                  u'ANN[*].CDNA_LEN',
                  u'ANN[*].CDS_POS',
                  u'ANN[*].CDS_LEN',
                  u'ANN[*].AA_POS',
                  u'ANN[*].AA_LEN',
                  u'ANN[*].DISTANCE',
                  u'ANN[*].ERRORS',
                  u'LOF[*].GENE',
                  u'LOF[*].GENEID',
                  u'LOF[*].NUMTR',
                  u'LOF[*].PERC',
                  u'NMD[*].GENE',
                  u'NMD[*].GENEID',
                  u'NMD[*].NUMTR',
                  u'NMD[*].PERC',
                  u'ANN[*].EFFECT',
                  u'ANN[*].IMPACT',
                  u'ANN[*].GENE',
                  u'ANN[*].GENEID',
                  u'ANN[*].FEATUREID',
                  u'ANN[*].BIOTYPE',
#                  u'dbNSFP_rs_dbSNP146',
                  u'dbNSFP_aapos',
                  u'dbNSFP_aaref',
                  u'dbNSFP_aaalt',
                  u'dbNSFP_Uniprot_acc_Polyphen2',
                  u'dbNSFP_Uniprot_id_Polyphen2',
                  u'dbNSFP_Uniprot_aapos_Polyphen2',
                  u'dbNSFP_1000Gp3_AF',
                  u'dbNSFP_ExAC_AF',
                  u'dbNSFP_Polyphen2_HVAR_pred',
                  u'dbNSFP_Polyphen2_HDIV_pred',
                  u'dbNSFP_CADD_phred',
                  u'dbNSFP_MetaSVM_pred',
                  u'dbNSFP_SIFT_pred',
                  u'dbNSFP_MetaLR_score',
                  u'dbNSFP_MetaLR_rankscore',
                  u'dbNSFP_MetaLR_pred',
                  u'dbNSFP_M_CAP_score',
                  u'dbNSFP_M_CAP_rankscore',
                  u'dbNSFP_M_CAP_pred',
                  u'syn_z',
                  u'mis_z',
                  u'lof_z',
                  u'pLI',
                  u'pRec',
                  u'pNull',
                  u'SFARIscore',
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
                  u'VARTYPE',
                  u'var_id',
                  u'v_id']


extra_cols = ['c_spark_genes',
              'c_cohort_freq',
              'coding_var',
              'c_effect_cat',
              'c_pop_freq',
              'c_allele_frac',
              'c_missense',
              'c_lof',
              'c_syn',
              'c_dmg_miss',
              'c_impact_lof']


def calcMetr(vn_df, msg=' '):
    tst = train.TrainTest('x',
                      '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt',
                      ['status'],
                      ['descr'])
    if vn_df.empty:
        return None
    vn_d = vn_df[~vn_df.var_id.duplicated()]
    tst.pred_y = numpy.array(vn_d['pred_labels'].astype(int))
    tst.test_set_y = tst.pred_y * 0
    tst.test_set_y[numpy.array(vn_d.status.isin(['Y']))] = 1
    print(msg)
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
                       prefix,
                       outp_dir,
                       config_file,
                       exac_anno='/mnt/xfs1/scratch/asalomatov/data/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt',
                       asd_gene_prob_anno='/mnt/xfs1/scratch/asalomatov/data/gene-scores/asd_gene_prediction_olga.csv',
                       ios_anno='/mnt/xfs1/scratch/asalomatov/data/gene-scores/ioss_lgd_rvis.scores.csv',
                       sfari_scores='/mnt/xfs1/scratch/asalomatov/data/SFARI/gene-score-only.csv'):
    with open(config_file, 'r') as f:
        cfg = yaml.safe_load(f)

    ped_file = cfg['ped_file']
    ped_file_extended = cfg['ped_file_extended']
    # populate ped DF
    if ped_file and ped_file_extended:
        sys.exit('only one of ped_file, ped_file_extended may be non-empty')
    if ped_file:
        myped = ped.Ped(ped_file)
        myped.addVcf(file_pat=cfg['vcf_pattern'])
        myped.ped.dropna(subset=['vcf'], inplace=True)
        myped.ped.reset_index(inplace=True)
    elif ped_file_extended:
        myped = ped.Ped(ped_file_extended, ['bam', 'vcf'])
    else:
        sys.exit('ped_file or ped_file_extended must be defined')
    exac = pandas.read_table(exac_anno)
#    asd_gene_prob_df = pandas.read_csv(asd_gene_prob_anno)
#    ios_anno_df = pandas.read_csv(ios_anno)
    sfari_scores_df = pandas.read_csv(sfari_scores)
    vn = pandas.read_table(infile)
    vn.columns = vn.columns.str.translate(None, '#')
    print(vn.shape)
    vn.ix[:, 'gene'] = vn['ANN[*].GENE']
    vn = vn.merge(
        exac[[u'syn_z', u'mis_z', u'lof_z', u'pLI', u'pRec', u'pNull', u'gene']],
        on='gene', how='left')
    vn = vn.merge(
        sfari_scores_df,
        on='gene', how='left')
    print(vn.shape)
    vn['v_id'] = vn.ind_id.astype(str) + '_' +\
                 vn['CHROM'].astype(str) + '_' +\
                 vn.POS.astype(str) + '_' +\
                 vn['ANN[*].GENE'] # + '_' +\
                 # vn['ANN[*].FEATUREID']
                 # vn['ANN[*].EFFECT'] + '_' +\
                 # vn['ANN[*].IMPACT']
    vn['var_id'] = vn.ind_id.astype(str) + '_' +\
                 vn['CHROM'].astype(str) + '_' +\
                 vn.POS.astype(str)
    vn['chr_pos'] = vn['CHROM'].astype(str) + '_' +\
                 vn.POS.astype(str)

    #vn = vn.merge(kv_vcf[['var_id', 'status']], on='var_id', how='left')
    print('before dedup')
    print(vn.shape)

    vn = vn[~vn.v_id.duplicated()]
    print('after dedup')
    print(vn.shape)
    
#    vn_all = vn
# stats before any filtering
#    print('\ndeduped and annotated vars, pred_labels value_counts:')
#    print(vn.pred_labels.value_counts())
#    print('deduped and annotated vars, test_labels value_counts:')
#    print(vn.status.value_counts())
#    calcMetr(vn, msg='deduped metrics')
#    vn_full = vn

    var_freq = vn.groupby('chr_pos').apply(lambda x: len(x['ind_id'].unique()))
    c_cohort_freq = var_freq > cfg['max_cohort_freq']
    var_freq_2 = var_freq[c_cohort_freq]
    # vn = vn[~vn.chr_pos.isin(var_freq_2.index)]
    vn['c_cohort_freq'] = ~vn.chr_pos.isin(var_freq_2.index)
#    print('\ncohort freq vars, pred_labels value_counts:')
#    print(vn.pred_labels.value_counts())
#    print('cohort freq vars, test_labels value_counts:')
#    print(vn.status.value_counts())
#    calcMetr(vn, msg='cohort_freq')
#    vn_diff =  getDiff(vn_full, vn, msg='cohort_freq')

    vn.ix[:, 'effect_cat'] = 'other'
    vn.ix[vn['ANN[*].EFFECT'].str.contains(
        '|'.join(cfg['snpeff']['effect_synon'])), 'effect_cat'] = 'syn' 
    vn.ix[vn['ANN[*].EFFECT'].str.contains(
        '|'.join(cfg['snpeff']['effect_dmgmis'])), 'effect_cat'] = 'mis' 
    vn.ix[vn['ANN[*].EFFECT'].str.contains(
        '|'.join(cfg['snpeff']['effect_lof'])), 'effect_cat'] = 'lof'
    vn['c_effect_cat'] = ~vn.effect_cat.isin(['other'])
#    print(vn.shape)
#    vn_full = vn
#    vn = vn.dropna(subset=['effect_cat'], axis=0)
#    print(vn.shape)

#    print('\neffects of interest vars, pred_labels value_counts:')
#    print(vn.pred_labels.value_counts())
#    print('effects of interest vars, test_labels value_counts:')
#    print(vn.status.value_counts())
#    calcMetr(vn, msg='effects metrics')
#    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='effects')])

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

#    vn_full = vn

    vn['c_pop_freq'] = (vn.dbNSFP_1000Gp3_AF < cfg['population_AF']) &\
                       (vn.dbNSFP_ExAC_AF < cfg['population_AF'])
#    vn = vn[(vn.dbNSFP_1000Gp3_AF < cfg['population_AF']) &
#            (vn.dbNSFP_ExAC_AF < cfg['population_AF'])]
#    print('\nAF vars, pred_labels value_counts:')
#    print(vn.pred_labels.value_counts())
#    print('AF vars, test_labels value_counts:')
#    print(vn.status.value_counts())
#    calcMetr(vn, msg='AF metrics')
#    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='pop_freq')])

#    vn_full = vn

    vn['c_biotype'] = vn['ANN[*].BIOTYPE'].str.contains(
        '|'.join(cfg['snpeff']['biotype']))
#    vn = vn[vn['ANN[*].BIOTYPE'].str.contains('|'.join(cfg['snpeff']['biotype']))]

#    print('\nprotein coding vars, pred_labels value_counts:')
#    print(vn.pred_labels.value_counts())
#    print('protein coding vars, test_labels value_counts:')
#    print(vn.status.value_counts())
#    calcMetr(vn, msg='protein coding metrics')
#    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='protein')])

#    vn_full = vn
    vn['allele_frac'] = vn.alt_DP.astype(float)/vn.DP
    vn['c_allele_frac'] = (vn.allele_frac > cfg['alt_allele_frac_range'][0]) &\
                          (vn.allele_frac < cfg['alt_allele_frac_range'][1])
#    vn = vn[(allele_frac > cfg['alt_allele_frac_range'][0]) & (allele_frac < cfg['alt_allele_frac_range'][1])]

#    print('\nallele fraction, pred_labels value_counts:')
#    print(vn.pred_labels.value_counts())
#    print('allel fraction vars, test_labels value_counts:')
#    print(vn.status.value_counts())
#    calcMetr(vn, msg='all fraction metrics')
#    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='allele_frac')])

    vn = vn.replace('ZZZ', '.')
#    if vn.empty:
#        print('No de novo mutation of interest')
#        return 0
    print('vn shape')
    print(vn.shape)
    vn['FILTER'] = vn.apply(func.getFieldFromVCF, args=(myped,), axis=1)
    vn = vn[~vn.FILTER.isnull()]
    print('vn shape in vcf')
    print(vn.shape)
    
    non_coding_vars = ['intron_variant', 'downstream_gene_variant',
                       'upstream_gene_variant', 'sequence_feature',
                       '5_prime_UTR_variant', '3_prime_UTR_variant']
    vn['coding_var'] = True
    for i in non_coding_vars:
        vn.ix[vn['ANN[*].EFFECT'] == i, 'coding_var'] = False

    c_missense = vn['effect_cat'] == 'mis'
    c_lof = vn['effect_cat'] == 'lof'
    c_syn = vn['effect_cat'] == 'syn'
    c_other = vn['effect_cat'] == 'other'
    vn['c_missense'] = c_missense
    vn['c_lof'] = c_lof
    vn['c_syn'] = c_syn
    
    print('M_CAP condition %s' % '|'.join(cfg['db_nsfp']['M_CAP_pred']))
    c_M_CAP_D = vn.dbNSFP_M_CAP_pred.str.contains(
        '|'.join(cfg['db_nsfp']['M_CAP_pred']))
    print(vn.dbNSFP_M_CAP_pred.value_counts())
    c_metaSVM_D = vn.dbNSFP_MetaSVM_pred.str.contains(
        '|'.join(cfg['db_nsfp']['metaSVM_pred']))
    print(vn.dbNSFP_CADD_phred.value_counts())
    c_cadd_null = vn.dbNSFP_CADD_phred.isin(['ZZZ', '.'])
    vn.ix[c_cadd_null, 'dbNSFP_CADD_phred'] = 0
    print(vn.dbNSFP_CADD_phred.value_counts())
    vn.ix[:, 'dbNSFP_CADD_phred'] = vn.dbNSFP_CADD_phred.astype(
        str).str.replace(',\.', ',0')
    print(vn.dbNSFP_CADD_phred.value_counts())
    vn.ix[:, 'dbNSFP_CADD_phred'] = vn.dbNSFP_CADD_phred.astype(
        str).str.replace('\.,', '0,')
    print(vn.dbNSFP_CADD_phred.value_counts())
    c_cadd_D = vn.dbNSFP_CADD_phred.astype(str).apply(
        lambda x: max(map(float, x.split(',')))) >=\
        cfg['db_nsfp']['cadd_phred']
    c_cadd_15 = vn.dbNSFP_CADD_phred.astype(str).apply(
        lambda x: max(map(float, x.split(',')))) >=\
        cfg['db_nsfp']['combined']['cadd_phred']
    c_poly_HVAR_D = vn.dbNSFP_Polyphen2_HVAR_pred.str.contains(
        '|'.join(cfg['db_nsfp']['combined']['polyphen2_pred']))
    c_poly_HDIV_D = vn.dbNSFP_Polyphen2_HVAR_pred.str.contains(
        '|'.join(cfg['db_nsfp']['combined']['polyphen2_pred']))
    c_sift_D = vn.dbNSFP_SIFT_pred.str.contains(
        '|'.join(cfg['db_nsfp']['combined']['sift_pred']))
    c_dmg_miss = c_M_CAP_D | c_metaSVM_D | c_cadd_D |\
                 ((c_poly_HDIV_D | c_poly_HVAR_D) & c_sift_D & c_cadd_15)
    print('N dmg_mis w/o M_CAP %s' % sum(c_metaSVM_D | c_cadd_D |
                                         ((c_poly_HDIV_D | c_poly_HVAR_D) &
                                          c_sift_D & c_cadd_15)))
    print('N dmg_mis w M_CAP %s' % sum(c_M_CAP_D | c_metaSVM_D | c_cadd_D |
                                       ((c_poly_HDIV_D | c_poly_HVAR_D) &
                                        c_sift_D & c_cadd_15)))
    vn['c_dmg_miss'] = c_dmg_miss
    c_impact_lof = vn['ANN[*].IMPACT'].str.contains(
        '|'.join(cfg['snpeff']['impact_lof']))
    vn['c_impact_lof'] = c_impact_lof
    c_all_denovo = vn.c_cohort_freq & vn.c_pop_freq &\
                   vn.c_allele_frac
    c_prev = vn.c_cohort_freq &\
             vn.c_pop_freq &\
             vn.c_allele_frac
             # vn.c_effect_cat &\
#             vn.c_biotype &\
    print('sum(c_prev)')
    print(sum(c_prev))
    c_spark_genes = vn['ANN[*].GENE'].str.contains(
        '|'.join(cfg['snpeff']['genes']))
    vn['c_spark_genes'] = c_spark_genes
    vn_mis = vn[c_dmg_miss & c_missense & c_prev]
    vn_mis_clinical = vn[c_dmg_miss & c_missense & c_prev & c_spark_genes]
    print('shape vn_mis')
    print(vn_mis.shape)
#    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn_mis, msg='dmg_miss')])
#    vn_full = vn[c_lof]
    vn_lof = vn[c_lof & c_impact_lof & c_prev]
    vn_lof_clinical = vn[c_lof & c_impact_lof & c_prev & c_spark_genes]
#    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn_lof, msg='impact_lof')])

    vn_syn = vn[c_syn & c_prev]
    vn_syn_clinical = vn[c_syn & c_prev & c_spark_genes]

    vn_other = vn[c_other & c_prev]
    vn_other_clinical = vn[c_other & c_prev & c_spark_genes]
#    print(vn.shape)

    c_FN = (vn.pred_labels == 0) & vn.status.isin(['Y'])
    vn_FN = vn[cols_to_output][c_FN]
    #vn_FN = vn_FN[~vn_FN.v_id.duplicated()]

    c_TP = (vn.pred_labels == 1) & vn.status.isin(['Y'])
    vn_TP = vn[cols_to_output][c_TP]
    #vn_TP = vn_TP[~vn_TP.v_id.duplicated()]
    
    var_type = cfg['variant_type']
    outp_suffix = '' # '{:%Y-%m-%d_%H-%M-%S-%f}'.format(datetime.datetime.now())
    
    def writeVariants(df, cols_to_output, var_type, prefix, suffix, outp_dir):
        if df.empty:
            print('%s is empty' % prefix)
            return None
        df = df[cols_to_output]
        df.columns = df.columns.str.replace('c_', '')
        df.to_csv(os.path.join(outp_dir,
                               '_'.join([prefix,
                                         var_type,
                                         suffix]) + '.csv'),
                  index=False)
    
    writeVariants(vn[c_all_denovo], cols_to_output[:-2] + extra_cols, var_type,
                  prefix, 'ALL_DENOVO', outp_dir)
    writeVariants(vn[vn.c_biotype], cols_to_output[:-2] + extra_cols, var_type,
                  prefix, 'ALL_DENOVO_CODING', outp_dir)
    writeVariants(vn_FN, cols_to_output[:-2] + extra_cols, var_type, prefix,
                  'FN', outp_dir)
    writeVariants(vn_TP, cols_to_output[:-2] + extra_cols, var_type, prefix,
                  'TP', outp_dir)
    writeVariants(vn_mis, cols_to_output[:-2] + extra_cols, var_type, prefix,
                  'MIS', outp_dir)
    writeVariants(vn_lof, cols_to_output[:-2] + extra_cols, var_type, prefix,
                  'LOF', outp_dir)
    writeVariants(vn_syn, cols_to_output[:-2] + extra_cols, var_type, prefix,
                  'SYN', outp_dir)
    writeVariants(vn_other, cols_to_output[:-2] + extra_cols, var_type, prefix,
                  'OTHER', outp_dir)
    writeVariants(vn_mis_clinical, cols_to_output[:-2] + extra_cols, var_type,
                  prefix + '_MIS', 'clinical', outp_dir)
    writeVariants(vn_lof_clinical, cols_to_output[:-2] + extra_cols, var_type,
                  prefix + '_LOF', 'clinical', outp_dir)
    writeVariants(vn_syn_clinical, cols_to_output[:-2] + extra_cols, var_type,
                  prefix + '_SYN', 'clinical', outp_dir)
    writeVariants(vn_other_clinical, cols_to_output[:-2] + extra_cols,
                  var_type, prefix + '_OTHER', 'clinical', outp_dir)


#    writeVariants(vn_diff, cols_to_output[:-2]+['step'], var_type,
#                  prefix + '_DIFF', outp_suffix, outp_dir)
        
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

    with open(os.path.join(outp_dir, 'cfg' + outp_suffix + '.yml'), 'w') as f:
        yaml.dump(cfg, f, default_flow_style=False)
    return vn[c_all_denovo]

if __name__ == '__main__':
    infile = sys.argv[1]
    pref = sys.argv[2]
    outp_dir = sys.argv[3]
    config_file = sys.argv[4]
    func.runInShell('mkdir -p ' + outp_dir)
    summarizeMutations(infile,
                       pref,
                       outp_dir,
                       config_file)



# summarizeMutations(infile,
#                       prefix,
#                       outp_dir,
#                       config_file,
#                      exac_anno='/mnt/scratch/asalomatov/data/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'):




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


