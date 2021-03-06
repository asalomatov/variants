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
import summarizeVariants

cols_to_output = summarizeVariants.cols_to_output[:] + ['HGVSc;Exon;Intron;HGVSp']
extra_cols = summarizeVariants.extra_cols[:]


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
                       infile_vep,
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


    #kv_vcf = pandas.read_csv('/mnt/xfs1/scratch/asalomatov/data/columbia/feature_sets/known/all_known.txt', sep='\t')
    #kv_vcf = kv_vcf[['ind_id','CHROM', 'POS', 'REF_offspring', 'ALT_base_offspring', 'status', 'descr', 'DP_offspring', 'DP_father', 'DP_mother']]
    #kv_vcf = kv_vcf[kv_vcf.descr.isin(['after'])]
    #kv_vcf['var_id'] = kv_vcf.ind_id.astype(str)+'_'+kv_vcf.CHROM.astype(str)+'_'+kv_vcf.POS.astype(str)
    #effects_of_interest = effects_loss_of_func + '|' + effect_damaging_missense + '|' + effect_synon
    exac = pandas.read_table(exac_anno)
    sfari_scores_df = pandas.read_csv(sfari_scores)
    vn = pandas.read_table(infile)
    vn.columns = vn.columns.str.translate(None, '#')
    vn.columns = [i.replace('[*]', '') for i in vn.columns]
    # read vep
    vep = func.readVcfToDF(infile_vep)
    vep = vep.merge(vep.apply(lambda row: func.vepVar2vcfVar(row, cfg['genome_ref']), axis=1),
                    right_index=True, left_index=True)

    vn.ix[:, 'gene'] = vn['ANN.GENE']
    vn = vn.merge(
        exac[[u'syn_z', u'syn_z_rank', u'syn_z_perc_rank',
              u'mis_z', u'mis_z_rank', u'mis_z_perc_rank',
              u'lof_z', u'lof_z_rank', u'lof_z_perc_rank',
              u'pLI', u'pLI_rank', u'pLI_perc_rank',
              u'pRec', u'pRec_rank', u'pRec_perc_rank',
              u'pNull', u'pNull_rank', u'pNull_perc_rank',
              u'gene']],
        on='gene', how='left')
    vn = vn.merge(
        sfari_scores_df,
        on='gene', how='left')
    print(vn.shape)
    vn['v_id'] = vn.ind_id.astype(str) + '_' +\
                 vn['CHROM'].astype(str) + '_' +\
                 vn.POS.astype(str) + '_' +\
                 vn['ANN.GENE'] # + '_' +\
                 # vn['ANN[*].FEATUREID']
                 # vn['ANN[*].EFFECT'] + '_' +\
                 # vn['ANN[*].IMPACT']
    vn['var_id'] = vn.ind_id.astype(str) + '_' +\
                 vn['CHROM'].astype(str) + '_' +\
                 vn.POS.astype(str)
    vn['chr_pos'] = vn['CHROM'].astype(str) + '_' +\
                 vn.POS.astype(str)
    vn['chr_pos_allel_tr'] = vn['CHROM'].astype(str) + '_' +\
                             vn.POS.astype(str) + '_' +\
                             vn.REF.astype(str) + '_' +\
                             vn.ALT.astype(str) + '_' +\
                             vn['ANN.FEATUREID'].astype(str)

    # if VEP anno is refseq based, translate to ensembl
    if vep.Feature.str.startswith('NM').sum() > 0 or\
                vep.Feature.str.startswith('XM').sum() > 0:
        ens_refseq = pandas.read_csv(cfg['ens_refseq'])
        vep['vep_transcript'] = vep.Feature.apply(
                lambda i: i.split('.')[0])
        vep = vep.merge(ens_refseq[['enst', 'nm']],
                        how='left',
                        left_on='vep_transcript',
                        right_on='nm')
    else:
        vep['enst'] = vep.Feature
    vep['chr_pos'] = vep['CHROM'].astype(str) + '_' +\
                             vep.POS.astype(str)
    vep['chr_pos_allel'] = vep['CHROM'].astype(str) + '_' +\
                             vep.POS.astype(str) + '_' +\
                             vep.REF.astype(str) + '_' +\
                             vep.ALT.astype(str)
    # join all transcripts
    vep_by_var = vep.groupby('chr_pos').apply(func.mergeFieldsForVariant).to_frame()
    vep_by_var.reset_index(inplace=True)
    vep_by_var.columns = ['chr_pos', 'HGVSc;Exon;Intron;HGVSp']
    vep['chr_pos_allel_tr'] = vep['CHROM'].astype(str) + '_' +\
                             vep.POS.astype(str) + '_' +\
                             vep.REF.astype(str) + '_' +\
                             vep.ALT.astype(str) + '_' +\
                             vep.enst.astype(str)
    print(vn.chr_pos_allel_tr)
    print(vep.chr_pos_allel_tr)
    print('vn dim before merging with vep:')
    print(vn.shape)
    vn = vn.merge(vep, how='left',
                  left_on='chr_pos_allel_tr',
                  right_on='chr_pos_allel_tr',
                  suffixes=['', '_vep'])
    vn = vn.merge(vep_by_var, how='left',
                  left_on='chr_pos',
                  right_on='chr_pos',
                  suffixes=['', '_vep'])
    # vn = vn.merge(kv_vcf[['var_id', 'status']], on='var_id', how='left')
    print('vn dim after merging with vep:')
    print(vn.shape)
    print('before dedup')
    print(vn.shape)
    vn_dups = vn[vn.v_id.duplicated()]
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
    vn.ix[vn['ANN.EFFECT'].str.contains(
        '|'.join(cfg['snpeff']['effect_synon'])), 'effect_cat'] = 'syn'
    vn.ix[vn['ANN.EFFECT'].str.contains(
        '|'.join(cfg['snpeff']['effect_dmgmis'])), 'effect_cat'] = 'mis'
    vn.ix[vn['ANN.EFFECT'].str.contains(
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

    vn['c_biotype'] = vn['ANN.BIOTYPE'].str.contains(
        '|'.join(cfg['snpeff']['biotype']))
#    vn = vn[vn['ANN.BIOTYPE'].str.contains('|'.join(cfg['snpeff']['biotype']))]

#    print('\nprotein coding vars, pred_labels value_counts:')
#    print(vn.pred_labels.value_counts())
#    print('protein coding vars, test_labels value_counts:')
#    print(vn.status.value_counts())
#    calcMetr(vn, msg='protein coding metrics')
#    vn_diff = pandas.concat([vn_diff, getDiff(vn_full, vn, msg='protein')])

#    vn_full = vn
    vn['allele_frac'] = None # vn.alt_DP.astype(float)/vn.DP
    vn['c_allele_frac'] = None # (vn.allele_frac > cfg['alt_allele_frac_range'][0]) &\
                        #  (vn.allele_frac < cfg['alt_allele_frac_range'][1])
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
    vn['FILTER'] = None # vn.apply(func.getFieldFromVCF, args=(myped,), axis=1)
#    vn = vn[~vn.FILTER.isnull()]
    print('vn shape in vcf')
    print(vn.shape)

    non_coding_vars = ['intron_variant', 'downstream_gene_variant',
                       'upstream_gene_variant', 'sequence_feature',
                       '5_prime_UTR_variant', '3_prime_UTR_variant']
    vn['coding_var'] = True
    for i in non_coding_vars:
        vn.ix[vn['ANN.EFFECT'] == i, 'coding_var'] = False

    c_missense = vn['effect_cat'] == 'mis'
    c_lof = vn['effect_cat'] == 'lof'
    c_syn = vn['effect_cat'] == 'syn'
    c_other = vn['effect_cat'] == 'other'
    vn['c_missense'] = c_missense
    vn['c_lof'] = c_lof
    vn['c_syn'] = c_syn

    c_M_CAP_D = vn.dbNSFP_M_CAP_pred.str.contains(
        '|'.join(cfg['db_nsfp']['M_CAP_pred']))
    c_metaSVM_D = vn.dbNSFP_MetaSVM_pred.str.contains(
        '|'.join(cfg['db_nsfp']['metaSVM_pred']))
#    c_metaSVM_null = vn.dbNSFP_MetaSVM_pred.isin(['ZZZ', '.'])

    c_cadd_null = vn.dbNSFP_CADD_phred.isin(['ZZZ', '.'])
    vn.ix[c_cadd_null, 'dbNSFP_CADD_phred'] = 0
    vn.ix[:, 'dbNSFP_CADD_phred'] = vn.dbNSFP_CADD_phred.astype(
        str).str.replace(',\.', ',0')
    vn.ix[:, 'dbNSFP_CADD_phred'] = vn.dbNSFP_CADD_phred.astype(
        str).str.replace('\.,', '0,')
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

#    c_sift_null = vn.dbNSFP_SIFT_pred.isin(['ZZZ', '.'])
    c_sift_D = vn.dbNSFP_SIFT_pred.str.contains(
        '|'.join(cfg['db_nsfp']['combined']['sift_pred']))
#   c_new = (vn.pred_labels == 1) & (~vn.status.isin(['Y']))

    c_dmg_miss = c_M_CAP_D | c_metaSVM_D | c_cadd_D |\
                 ((c_poly_HDIV_D | c_poly_HVAR_D) & c_sift_D & c_cadd_15)
    print('N dmg_mis w/o M_CAP %s' % sum(c_metaSVM_D | c_cadd_D |
                                         ((c_poly_HDIV_D | c_poly_HVAR_D) &
                                          c_sift_D & c_cadd_15)))
    print('N dmg_mis w M_CAP %s' % sum(c_M_CAP_D | c_metaSVM_D | c_cadd_D |
                                       ((c_poly_HDIV_D | c_poly_HVAR_D) &
                                        c_sift_D & c_cadd_15)))
    vn['c_dmg_miss'] = c_dmg_miss
    vn['c_dmg_miss_woMCAP'] = c_metaSVM_D | c_cadd_D |\
                              ((c_poly_HDIV_D | c_poly_HVAR_D) &
                               c_sift_D & c_cadd_15)

#   vn_full = vn[c_missense]
    c_impact_lof = vn['ANN.IMPACT'].str.contains(
        '|'.join(cfg['snpeff']['impact_lof']))
    vn['c_impact_lof'] = c_impact_lof

    c_all_denovo = vn.c_cohort_freq & vn.c_pop_freq
    c_prev = vn.c_cohort_freq &\
             vn.c_pop_freq
#             vn.c_allele_frac
             # vn.c_effect_cat &\
#             vn.c_biotype &\
    print('sum(c_prev)')
    print(sum(c_prev))
    c_spark_genes = vn['ANN.GENE'].str.contains(
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
    vn_lof_dmis_clinical = pandas.concat([vn_lof_clinical, vn_mis_clinical])

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
        df.columns = map(lambda i: i[2:] if i.startswith('c_') else i,
                         df.columns)
        df.to_csv(os.path.join(outp_dir,
                               '_'.join([prefix,
                                         var_type,
                                         suffix]) + '.csv'),
                  index=False)

    writeVariants(vn[c_all_denovo], cols_to_output + extra_cols, var_type,
                  prefix, 'ALL_DENOVO', outp_dir)
    writeVariants(vn[vn.c_biotype], cols_to_output + extra_cols, var_type,
                  prefix, 'ALL_DENOVO_CODING', outp_dir)
    writeVariants(vn_FN, cols_to_output + extra_cols, var_type, prefix,
                  'FN', outp_dir)
    writeVariants(vn_TP, cols_to_output + extra_cols, var_type, prefix,
                  'TP', outp_dir)
    writeVariants(vn_mis, cols_to_output + extra_cols, var_type, prefix,
                  'MIS', outp_dir)
    writeVariants(vn_lof, cols_to_output + extra_cols, var_type, prefix,
                  'LOF', outp_dir)
    writeVariants(vn_syn, cols_to_output + extra_cols, var_type, prefix,
                  'SYN', outp_dir)
    writeVariants(vn_other, cols_to_output + extra_cols, var_type, prefix,
                  'OTHER', outp_dir)
    writeVariants(vn_mis_clinical, cols_to_output + extra_cols, var_type,
                  prefix + '_MIS', 'clinical', outp_dir)
    writeVariants(vn_lof_clinical, cols_to_output + extra_cols, var_type,
                  prefix + '_LOF', 'clinical', outp_dir)
    writeVariants(vn_lof_dmis_clinical, cols_to_output + extra_cols, var_type,
                  prefix + '_LOF_DMIS', 'clinical', outp_dir)
    writeVariants(vn_syn_clinical, cols_to_output + extra_cols, var_type,
                  prefix + '_SYN', 'clinical', outp_dir)
    writeVariants(vn_other_clinical, cols_to_output + extra_cols,
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
#                      exac_anno='/mnt/xfs1/scratch/asalomatov/data/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'):




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
    #       'vcf_pattern': '/mnt/xfs1/scratch/asalomatov/data/columbia/vcf/%s_%s-02_%s-01.annotated-norm.vcf.gz',
    #       'bam_pattern': '/mnt/xfs1/scratch/asalomatov/data/columbia/bam/Sample.%s.bam',
    #       'ped_file': '/mnt/xfs1/scratch/asalomatov/data/columbia/pcgc_ped.txt',
    #       'known_variants': '/mnt/xfs1/scratch/asalomatov/data/columbia/pcgc_denovo_snp.tsv',
    #       'output_directory': '/mnt/xfs1/scratch/asalomatov/data/columbia/feature_sets_01'
    #}


