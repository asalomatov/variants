import sys, os
import pandas
from variants import func
import ped
import argparse
import collections
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants/work/')
from init_baylor_info import spID2labID
from init_baylor_info import labID2spID
from init_baylor_info import lab2batch


sf_dn_c_columns = ['c_cohort_freq', 'c_biotype', 'c_effect_cat',
                   'c_pop_freq', 'c_allele_frac', 'c_missense', 'c_lof',
                   'c_syn', 'c_dmg_miss', 'c_impact_lof']
# c_cohort_freq,c_biotype,c_effect_cat,c_pop_freq,c_allele_frac,c_missense,c_lof,c_syn,c_dmg_miss,c_impact_lof


cols_to_output = [u'SP_id',
                  u'lab_id',
                  u'CHROM',
                  u'POS',
                  u'ID',
                  u'REF',
                  u'ALT',
                  u'ANN[*].GENE',
                  u'SFARIscore',
                  u'ANN[*].EFFECT',
                  u'ANN[*].IMPACT',
                  u'in_hc',
                  u'in_fb',
                  u'in_pl',
                  u'SF', u'Reg', u'Clmb',
                  u'pred_prob',
                  u'ref_DP',
                  u'alt_DP',
                  u'DP',
                  u'allele_frac',
                  u'DP_offspring',
                  u'DP_father',
                  u'DP_mother',
                  u'ANN[*].GENEID',
                  u'ANN[*].FEATUREID',
                  u'ANN[*].BIOTYPE',
                  u'effect_cat',
                  u'dbNSFP_rs_dbSNP146',
                  u'dbNSFP_aapos',
                  u'dbNSFP_aaref',
                  u'dbNSFP_aaalt',
                  u'dbNSFP_Uniprot_acPolyphen2',
                  u'dbNSFP_Uniprot_id_Polyphen2',
                  u'dbNSFP_Uniprot_aapos_Polyphen2',
                  u'dbNSFP_1000Gp3_AF',
                  u'dbNSFP_ExAC_AF',
                  u'dbNSFP_Polyphen2_HVAR_pred',
                  u'dbNSFP_Polyphen2_HDIV_pred',
                  u'dbNSFP_CADD_phred',
                  u'dbNSFP_MetaSVM_pred',
                  u'dbNSFP_SIFT_pred',
                  u'syn_z', u'mis_z',
                  u'lof_z', u'pLI', u'pRec', u'pNull',
                  u'spidex_dpsi_max_tissue',
                  u'spidex_dpsi_zscore', u'spidex_gene', u'spidex_strand',
                  u'spidex_transcript', u'spidex_exon_number',
                  u'spidex_location',
                  u'spidex_cds_type', u'spidex_ss_dist', u'FILTER',
                  u'spark_genes',
                  u'asd_gene_score',
                  u'LDGscore', u'LGDrank',
                  u'RVIS', u'RVISrank', u'RVIS_LGDrank',
                  u'cohort_freq', u'coding_var', u'effect_cat.1', u'pop_freq',
                  u'allele_frac.1', u'missense', u'lof', u'syn', u'dmg_miss',
                  u'impact_lof', u'var_type', u'status',  u'batch', u'v_id']


def varType(row):
    c1 = len(row['REF']) == len(row['ALT'])
    if c1:
        return 'SNP'
    else:
        return 'INDEL'


def isPossibleDeNovo(row, clr):
    if len(row['REF']) == len(row['ALT']):
        v_t = 'snp'
    else:
        v_t = 'indels'
    cmd = ' '.join(['cat',
                    os.path.join(
                        sf_calls_dir,
                        v_t,
                        clr,
                        row['ind_id']),
                    '|',
                    'grep',
                    row['CHROM'],
                    '|',
                    'grep',
                    str(row['POS'])])
    res = func.runInShell(cmd)
    return not bool(res)


def isDeNovo(row, clr):
    if len(row['REF']) == len(row['ALT']):
        v_t = 'snp'
        v_tt = 'SNP'
    else:
        v_t = 'indels'
        v_tt = 'INDEL'
    cmd = ' '.join(['cat',
                    os.path.join(
                        sf_calls_dir,
                        v_t,
                        clr,
                        row['batch'],
                        row['ind_id'] + '-' + v_tt +
                        '-class.csv'),
                    '|',
                    'grep',
                    row['ind_id'],
                    '|',
                    'grep',
                    row['CHROM'],
                    '|',
                    'grep',
                    str(row['POS'])])
    res = func.runInShell(cmd)
    return not bool(res)


def whyNotDeNovo(row, clr):
    cmd = ' '.join(['cat',
                    os.path.join(sf_calls_dir, clr, clr+'-SNP-class_ALL_SNP__.csv'),
                    '|',
                    'grep',
                    row['ind_id'],
                    '|',
                    'grep',
                    row['CHROM'],
                    '|',
                    'grep',
                    str(row['POS'])])
    res = func.runInShell(cmd, True)
    print res
    if row['is_dn'] and not row['is_filt']:
        res = res.split('\n')[0]
        outp = res.split(',')[-len(sf_dn_c_columns):]
    else:
        outp = [None] * len(sf_dn_c_columns)

    ser_out = pandas.Series(outp, sf_dn_c_columns, dtype=str)
    print ser_out
    # ser_out = ser_out[sf_dn_c_columns]
    return ser_out


def isFiltered(row, clr):
    cmd = ' '.join(['cat',
                    os.path.join(sf_calls_dir, clr, clr+'-SNP-class_MIS_SNP__.csv'),
                    os.path.join(sf_calls_dir, clr, clr+'-SNP-class_SYN_SNP__.csv'),
                    os.path.join(sf_calls_dir, clr, clr+'-SNP-class_LOF_SNP__.csv'),
                    '|',
                    'grep',
                    row['ind_id'],
                    '|',
                    'grep',
                    row['CHROM'],
                    '|',
                    'grep',
                    str(row['POS'])])
    res = func.runInShell(cmd)
    return not bool(res)


def isInVcf(row, ped_obj):
    x = func.getFieldFromVCF(row, ped_obj)
    if x is None:
        return False
    else:
        return True


def isInVcfSNP(df):
    res = df.in_vcf_hc |\
          df.in_vcf_fb |\
          df.in_vcf_pl
    return res


def isInVcfINDEL(df):
    res = df.in_vcf_HC &\
          df.in_vcf_JHC &\
          df.in_vcf_FB &\
          df.in_vcf_PL
    return res


def isPossibleDeNovoSNP(df):
    res = df.possib_dn_HC |\
          df.possib_dn_JHC |\
          df.possib_dn_FB |\
          df.possib_dn_PL
    return res


def isPossibleDeNovoINDEL(df):
    res = df.possib_dn_HC &\
          df.possib_dn_JHC &\
          df.possib_dn_FB &\
          df.possib_dn_PL
    return res


def isDeNovoSNP(df):
    res = df.is_dn_hc |\
          df.is_dn_fb |\
          df.is_dn_pl
    return res


def isDeNovoINDEL(df):
    res = df.is_dn_HC &\
          df.is_dn_JHC &\
          df.is_dn_FB &\
          df.is_dn_PL
    return res


def isFilteredSNP(df):
    res = df.is_filt_HC |\
          df.is_filt_JHC |\
          df.is_filt_FB |\
          df.is_filt_PL
    return res


def get_spID(x, lab2sp_dict):
    if x[:2] == 'SP':
        return x
    else:
        return lab2sp_dict[x]


if __name__ == '__main__':
    print sys.argv
    arg_parser = argparse.ArgumentParser(
        description='Create a batch script for igv,\
        run igv in batch mode to generate\
        snapshots')
    arg_parser.add_argument('variant_callers',
                            help='a comma delimeted string\
                            such as hc,fb,pl')
    arg_parser.add_argument('batch_dir',
                            help='a coma delimited string, eg b1,b3,b4')
    arg_parser.add_argument('--reg_dn_calls',
                            default='/mnt/ceph/users/asalomatov/regeneron_spark_pilot/de_novo_variants/spark_regeneron_pilot_24_de_novo_calls.csv',
                            help='has ind_id, CHROM, POS, REF, ALT, gene')
    arg_parser.add_argument('--sy_dn_calls',
                            default='none',
                            help='lab_id, CHROM, POS, REF, ALT, DESCR,  GENE')
    arg_parser.add_argument('--data_dir',
                            default='/mnt/scratch/asalomatov/regeneron_spark_pilot/denovo/',
                            type=str,
                            help='directory where SF data resides')
    arg_parser.add_argument('--seq_center',
                            default='reg',
                            type=str,
                            help='bay, or reg')
    arg_parser.add_argument('--gen_other',
                            default='Yes',
                            type=str,
                            help='Yes, or No')
    arg_parser.add_argument('--extract_batch',
                            default='b6',
                            type=str,
                            help='batch to output along with the full set')
    args = arg_parser.parse_args()
    print args
    asd_gene_prob_anno = '/mnt/scratch/asalomatov/data/gene-scores/asd_gene_prediction_olga.csv'
    ios_anno = '/mnt/scratch/asalomatov/data/gene-scores/ioss_lgd_rvis.scores.csv'
    asd_gene_prob_df = pandas.read_csv(asd_gene_prob_anno)
    asd_gene_prob_df = asd_gene_prob_df[['symbol', 'score']]
    asd_gene_prob_df.columns = ['gene', 'asd_gene_score']
    ios_anno_df = pandas.read_csv(ios_anno)
    ios_anno_df.columns = ['gene', 'LDGscore', 'LGDrank',  'RVIS', 'RVISrank', 'RVIS_LGDrank']
    callers = args.variant_callers.split(',')
    sf_calls_dir = args.data_dir
    sf_seq_center = args.seq_center
    batch_dirs = args.batch_dir.split(',')
    calls_list = []
    for batch_dir in batch_dirs:
        # read data sets
        sf_calls_dict = collections.OrderedDict()
        clr_dict = {}
        for clr in callers:
            snps = pandas.read_csv(
                os.path.join(sf_calls_dir,
                             'snp',
                             clr,
                             batch_dir,
                             '%s-SNP-class_SNP_ALL_DENOVO.csv' % clr),
                dtype=str)
            snps['var_type'] = 'SNP'
            indels = pandas.read_csv(
                os.path.join(sf_calls_dir,
                             'indels',
                             clr,
                             batch_dir,
                             '%s-INDEL-class_INDEL_ALL_DENOVO.csv' % clr),
                dtype=str)
            indels['var_type'] = 'INDEL'
            df = pandas.concat([snps, indels])
            df.reset_index(inplace=True, drop=True)
            df['var_id'] = df.ind_id.astype(str) + '_' +\
                           df.CHROM.astype(str) + '_' +\
                           df.POS.astype(str) + '_' +\
                           df.var_type
            clr_dict[clr] = set(df.var_id)
            sf_calls_dict[clr] = df
        sf_calls = pandas.concat(sf_calls_dict)
        sf_calls = sf_calls[~sf_calls.var_id.duplicated()]
        sf_calls.reset_index(drop=True, inplace=True)
        for clr in callers:
            sf_calls['in_' + clr] = sf_calls.var_id.apply(lambda x:
                                                          x in clr_dict[clr])
        sf_calls['SP_id'] = sf_calls.ind_id  # .apply(lambda i:
        # get_spID(i, labID2spID))
        sf_calls['lab_id'] = sf_calls.SP_id  # .apply(lambda z: spID2labID[z])
        sf_calls['batch'] = batch_dir
        calls_list.append(sf_calls)

    sf_calls = pandas.concat(calls_list)
    sf_calls['chr_pos_type'] = sf_calls.CHROM.astype(str) + '_' +\
                           sf_calls.POS.astype(str) + '_' +\
                           sf_calls.var_type

#    sf_calls['lab_id'] = sf_calls.ind_id

    dpls = list(sf_calls.chr_pos_type[sf_calls.chr_pos_type.duplicated()])
    print 'chr_pos dups eliminated: ', len(dpls)
    sf_calls = sf_calls[~sf_calls.chr_pos_type.isin(dpls)]
    sf_calls['v_id'] = sf_calls.lab_id.astype(str) + '_' +\
                       sf_calls.CHROM.astype(str) + '_' +\
                       sf_calls.POS.astype(str)
    if args.reg_dn_calls != 'none':
        reg_calls = pandas.read_csv(args.reg_dn_calls, dtype='str')
        reg_calls['lab_id'] = reg_calls.ind_id.apply(lambda i: i.split('_')[1])
        reg_calls['SP_id'] = reg_calls.lab_id
        c_indel = reg_calls.ALT == '-'  #  reg_calls.POS != reg_calls.END
        reg_calls.POS[c_indel] = reg_calls.POS.astype(int) - 1
        reg_ref_atpos = reg_calls.apply(lambda row: func.refAtPos(row['CHROM'], row['POS'], '/mnt/scratch/asalomatov/data/b38/genome.fa'), axis=1)
        reg_calls.ALT[c_indel] = reg_ref_atpos[c_indel]
        reg_calls.REF[c_indel] = reg_calls.ALT[c_indel] + reg_calls.REF[c_indel]
        reg_calls['v_id'] = reg_calls.lab_id.astype(str) + '_' +\
                            reg_calls.CHROM.astype(str) + '_' +\
                            reg_calls.POS.astype(str)
        reg_set = set(reg_calls.v_id)
    else:
        reg_set = set()
    if args.sy_dn_calls != 'none':
        sy_calls = pandas.read_csv(args.sy_dn_calls, dtype='str')
        sy_calls['v_id'] = sy_calls.lab_id.astype(str) + '_' +\
                           sy_calls.CHROM.astype(str) + '_' +\
                           sy_calls.POS.astype(str)
        sy_set = set(sy_calls.v_id)
    else:
        sy_set = set()
    sf_set = set(sf_calls.v_id)
    sf_set_nonsyn = set(sf_calls.v_id[sf_calls.effect_cat != 'syn'])
    sf_calls['SF'] = True
#    sf_calls['SP_id'] = sf_calls.ind_id
    sf_calls['Clmb'] = sf_calls.v_id.apply(lambda z: z in sy_set)
    sf_calls['Reg'] = sf_calls.v_id.apply(lambda z: z in reg_set)
    sf_calls['sort_cat'] = 10
    sf_calls.ix[sf_calls.impact_lof == 'True', 'sort_cat'] = 1
    sf_calls.ix[sf_calls.dmg_miss == 'True', 'sort_cat'] = 2
    sf_calls.ix[(sf_calls.coding_var == 'True') & (sf_calls.dmg_miss == 'False') &
                (sf_calls.missense == 'True'), 'sort_cat'] = 3
    sf_calls.ix[(sf_calls.coding_var == 'True') &  (sf_calls.sort_cat > 3),
                'sort_cat'] = 4
    sf_calls_research = sf_calls.sort_values(['sort_cat', 'batch', 'lab_id',
                                              'CHROM', 'POS'])
    print('sf_calls_research dim:')
    print(sf_calls_research.shape)
    sf_calls_research['gene'] = sf_calls_research['ANN[*].GENE']
    sf_calls_research = sf_calls_research.merge(asd_gene_prob_df,
                                                on='gene', how='left')
    print('sf_calls_research dim:')
    print(sf_calls_research.shape)
    sf_calls_research = sf_calls_research.merge(ios_anno_df,
                                                on='gene', how='left')
    print('sf_calls_research dim:')
    print(sf_calls_research.shape)
    c_pass_hard_filters = sf_calls_research.FILTER.isin(['.', 'PASS'])
    sf_calls_research.ix[~c_pass_hard_filters, cols_to_output + ['sort_cat']].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_all_research_failed_hfilt.csv',
        index=False)
    sf_calls_research = sf_calls_research.ix[c_pass_hard_filters, :]
    sf_calls_research.reset_index(drop=True, inplace=True)
    sf_calls_research_failed = sf_calls_research.ix[~c_pass_hard_filters, :]
    sf_calls_research_failed.reset_index(drop=True, inplace=True)
    sf_calls_research[cols_to_output + ['sort_cat']].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_all_research.csv',
        index=False)
    sf_dmiss = (sf_calls_research.missense.astype(str) == 'True') & (sf_calls_research.dmg_miss.astype(str) == 'True')
    sf_lof = (sf_calls_research.lof.astype(str) == 'True') & (sf_calls_research.impact_lof.astype(str) == 'True')
    sf_spark = (sf_calls_research.spark_genes.astype(str) == 'True')
    sf_calls_research[cols_to_output + ['sort_cat']][sf_lof].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_lof.csv',
        index=False)
    sf_calls_research[cols_to_output + ['sort_cat']][sf_dmiss].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_dmis.csv',
        index=False)
    sf_calls_research[cols_to_output + ['sort_cat']][sf_lof].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_lof.csv',
        index=False)
    sf_calls_research[cols_to_output + ['sort_cat']][sf_spark].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_clinical.csv',
        index=False)
    sf_calls_research[cols_to_output + ['sort_cat']][sf_spark & (sf_lof | sf_dmiss)].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_clinical_lofdmis.csv',
        index=False)
    sf_calls_research[cols_to_output + ['sort_cat']][sf_calls_research.missense.astype(str) == 'True'].to_csv(
        sf_calls_dir + '/SF_denovo_' + '_'.join(batch_dirs) + '_missense.csv',
        index=False)

#    sys.exit('stop!')


    # now all other calls 
    if args.gen_other == 'Yes':
        other_calls = reg_calls
        other_calls.ind_id = other_calls.SP_id
        #pandas.concat([reg_calls,
        #                             bcm_calls[codif_calls.columns],
        #                             sy_calls[codif_calls.columns]])
        other_calls = other_calls[~other_calls.v_id.duplicated()]
        other_calls['in_SF'] = other_calls.v_id.apply(lambda i: i in sf_set)
        other_calls = other_calls[~other_calls.in_SF]
        other_calls['in_Clmb'] = other_calls.v_id.apply(lambda i:
                                                        i in sy_set)
        other_calls['batch'] = 'b1' #  'b' + other_calls.lab_id.apply(lambda i:
                                    #                      str(lab2batch[i]))
        #other_calls.ix[other_calls.batch.isin(['b1', 'b2']), 'batch'] = 'b1-2'
        #other_calls['ind_id'] = other_calls.lab_id.astype(str)
        #other_calls.ix[other_calls.batch != 'b1-2', 'ind_id'] =\
        #    other_calls.lab_id[other_calls.batch != 'b1-2'].apply(
        #        lambda i: get_spID(i, labID2spID))
        for clr in callers:
            ped_by_btch = {}
            for batch_dir in batch_dirs:
                if sf_seq_center == 'bay':
                    if batch_dir == 'b1-2':
                        path_toPed = '/mnt/scratch/asalomatov/data/SPARK/ped/spark_%s.ped' % clr
                    else:
                        path_toPed = '/mnt/scratch/asalomatov/data/SPARK/ped/spark_spID_%s_ext_%s.ped' % (batch_dir, clr)
                elif sf_seq_center == 'reg':
                    path_toPed = '/mnt/ceph/users/asalomatov/regeneron_spark_pilot/ped/spark_%s_ext_%s.ped' % (batch_dir, clr)  
                else:
                    sys.exit('unknown seq_center')
                myped = ped.Ped(path_toPed, ['bam', 'vcf'])
                ped_by_btch[batch_dir] = myped
                # check if variants are present in vcf files
            other_calls['in_vcf_' + clr] = other_calls.apply(
                lambda row:
                isInVcf(row,
                        ped_by_btch[row['batch']]),
                axis=1)
            other_calls['is_dn_' + clr] = other_calls.apply(
                lambda row:
                isDeNovo(row, clr),
                axis=1)
        
        other_calls['inVCF'] = isInVcfSNP(other_calls)
        other_calls['isDeNovo'] = isDeNovoSNP(other_calls)
        other_calls.to_csv(os.path.join(sf_calls_dir,
                                        'other_calls_' + '_'.join(batch_dirs) + '.csv'),
                           index=False)
    
        sys.exit('done, annotate other calls')

    other_ann = pandas.read_csv(
        os.path.join(
            sf_calls_dir,
            'other_calls_' + '_'.join(batch_dirs)) +
        '_ALL_ALL_DENOVO.csv')
#    other_ann['SP_id'] = other_ann.ind_id.apply(lambda z: get_spID(z, labID2spID))
#    other_ann['lab_id'] = other_ann.SP_id.apply(lambda z: spID2labID[z])
    other_ann['v_id'] = other_ann.ind_id.astype(str) + '_' +\
                        other_ann.CHROM.astype(str) + '_' +\
                        other_ann.POS.astype(str)
    other_ann['var_type'] = other_ann.apply(varType, axis=1)
    other_ann['var_id'] = other_ann.ind_id.astype(str) + '_' +\
                              other_ann.CHROM.astype(str) + '_' +\
                              other_ann.POS.astype(str) + '_' +\
                              other_ann.var_type
    other_ann['Reg'] = other_ann.v_id.apply(lambda z: z in reg_set)
    other_ann['Clmb'] = other_ann.v_id.apply(lambda z: z in sy_set)
    other_ann['SF'] = other_ann.v_id.apply(lambda z: z in sf_set)
    other_ann['batch'] = 'b1'  #  + other_ann.lab_id.apply(lambda i: str(lab2batch[i]))
    other_ann['sort_cat'] = 10
    other_ann.ix[other_ann.impact_lof, 'sort_cat'] = 1
    other_ann.ix[other_ann.dmg_miss, 'sort_cat'] = 2
    other_ann.ix[(other_ann.coding_var) & (other_ann.dmg_miss) &
                 (other_ann.missense), 'sort_cat'] = 3
    other_ann.ix[(other_ann.coding_var) & (other_ann.sort_cat > 3),
                 'sort_cat'] = 4
    other_ann['gene'] = other_ann['ANN[*].GENE']
    other_ann = other_ann.merge(asd_gene_prob_df, on='gene', how='left')
    other_ann = other_ann.merge(ios_anno_df, on='gene', how='left')

    all_calls_research = pandas.concat([sf_calls_research, other_ann])
    all_calls_research = all_calls_research.sort_values(['sort_cat', 'batch', 'lab_id',
                                            'CHROM', 'POS'])
    all_calls_research['N_centers'] = all_calls_research.SF.astype(int) +\
                                      all_calls_research.Reg.astype(int) +\
                                      all_calls_research.Clmb.astype(int)
    all_calls_research['var_type'] = all_calls_research.apply(varType, axis=1)
    all_calls_research.SP_id = all_calls_research.ind_id
    all_calls_research.lab_id = all_calls_research.ind_id
    def concatCenters(row):
        res = []
        if row['SF']: res.append('SF')
        if row['Reg']: res.append('Reg')
        if row['Clmb']: res.append('Clmb')
        return '_'.join(res)
    all_calls_research['centers'] = all_calls_research.apply(
        concatCenters, axis=1)
    all_calls_research.columns = [i.replace('ANN[*].', '') for i in all_calls_research.columns]
    cols_to_output = cols_to_output + ['centers', 'N_centers']
    cols_to_output = [i.replace('ANN[*].', '') for i in cols_to_output]
    all_lof = (all_calls_research.lof.astype(str) == 'True') & (all_calls_research.impact_lof.astype(str) == 'True')
    all_dmis = (all_calls_research.missense.astype(str) == 'True') & (all_calls_research.dmg_miss.astype(str) == 'True')
    all_mis = (all_calls_research.missense.astype(str) == 'True')
    all_spark = (all_calls_research.spark_genes.astype(str) == 'True') 
    all_lofdmiss = all_calls_research[all_lof | all_dmis]
    all_calls_research[cols_to_output + ['sort_cat']].to_csv(
        os.path.join(sf_calls_dir,
                     'all_calls_' + '_'.join(batch_dirs) + '_all_research.csv'),
        index=False)
    all_calls_research[cols_to_output + ['sort_cat']][all_lof].to_csv(
        os.path.join(sf_calls_dir,
                     'all_calls_' + '_'.join(batch_dirs) + '_lof.csv'),
        index=False)
    all_calls_research[cols_to_output + ['sort_cat']][all_dmis].to_csv(
        os.path.join(sf_calls_dir,
                     'all_calls_' + '_'.join(batch_dirs) + '_dmis.csv'),
        index=False)
    all_calls_research[cols_to_output + ['sort_cat']][all_mis].to_csv(
        os.path.join(sf_calls_dir,
                     'all_calls_' + '_'.join(batch_dirs) + '_mis.csv'),
        index=False)
    all_calls_research[cols_to_output + ['sort_cat']][all_spark].to_csv(
        os.path.join(sf_calls_dir,
                     'all_calls_' + '_'.join(batch_dirs) + '_clinical.csv'),
        index=False)
    all_calls_research[cols_to_output + ['sort_cat']][
        all_spark & (all_lof | all_dmis)].to_csv(
        os.path.join(sf_calls_dir,
                     'all_calls_' + '_'.join(batch_dirs) + '_clinical_lofdmis.csv'),
        index=False)
    c_b = args.extract_batch == all_calls_research.batch
    if sum(c_b) == 0:
        print 'no single batch to output'
    else:
        all_calls_research[cols_to_output + ['sort_cat']][c_b].to_csv(
            os.path.join(sf_calls_dir,
                         'all_calls_' + args.extract_batch + '_all_research.csv'),
            index=False)
        all_calls_research[cols_to_output + ['sort_cat']][all_lof & c_b].to_csv(
            os.path.join(sf_calls_dir,
                         'all_calls_' + args.extract_batch + '_lof.csv'),
            index=False)
        all_calls_research[cols_to_output + ['sort_cat']][all_dmis & c_b].to_csv(
            os.path.join(sf_calls_dir,
                         'all_calls_' + args.extract_batch + '_dmis.csv'),
            index=False)
        all_calls_research[cols_to_output + ['sort_cat']][all_mis & c_b].to_csv(
            os.path.join(sf_calls_dir,
                         'all_calls_' + args.extract_batch + '_mis.csv'),
            index=False)
        all_calls_research[cols_to_output + ['sort_cat']][all_spark & c_b].to_csv(
            os.path.join(sf_calls_dir,
                         'all_calls_' + args.extract_batch + '_clinical.csv'),
            index=False)
        all_calls_research[cols_to_output + ['sort_cat']][
            all_spark & (all_lof | all_dmis) & c_b].to_csv(
                os.path.join(sf_calls_dir,
                             'all_calls_' + args.extract_batch + '_clinical_lofdmis.csv'),
                index=False)



