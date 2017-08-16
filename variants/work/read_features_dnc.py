import pandas


def concatClr(row):
    res = []
    if row['fb']: res.append('fb')
    if row['hc']: res.append('hc')
    return '_'.join(res)


def read_candidates(clr, vartype):
    ssc_cand = pandas.read_table('/mnt/scratch/asalomatov/data/SSC/wes/de-novo-class/%(clr)s/%(vartype)s/all_%(vartype)s_features.tsv' % locals()) 
    ssc_known = pandas.read_table('/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_denovo_clean_%(vartype)s.tsv' % locals())
    ssc_known = ssc_known[~ssc_known.duplicated()]
    ssc_known['chr_pos'] = ssc_known.chr.astype(str) + '_' + ssc_known.pos.astype(str)
    ssc_known['v_id'] = ssc_known.ind_id.astype(str) + '_' + ssc_known.chr.astype(str) + '_' + ssc_known.pos.astype(str)

#    ssc_known = ssc_known[ssc_known.descr != 'Iossifov']
    print ssc_known.descr.value_counts()
    print ssc_known.status.value_counts()
    ssc_cand.shape
    ssc_cand['v_id'] = ssc_cand.ind_id.astype(str) + '_' + ssc_cand.CHROM.astype(str) + '_' + ssc_cand.POS.astype(str)
    ssc_cand['chr_pos'] = ssc_cand.CHROM.astype(str) + '_' + ssc_cand.POS.astype(str)
    ssc_cand_freq_2 = ssc_cand.chr_pos.value_counts()
    snp_freq = ssc_cand_freq_2.index[ssc_cand_freq_2 > 2].tolist()
    print '> 2', len(snp_freq)
    snp_freq = ssc_cand_freq_2.index[ssc_cand_freq_2 > 1].tolist()
    print '> 1', len(snp_freq)
    print 'before chr_pos dedup', ssc_cand.shape
    ssc_cand = ssc_cand[~ssc_cand.chr_pos.isin(snp_freq)]
    print 'after chr_pos dedup', ssc_cand.shape
    ssc_cand = ssc_cand.merge(ssc_known[['descr', 'status', 'v_id']], how='left',
                              on='v_id')
    ssc_cand['vartype'] = vartype
    return ssc_cand


if __name__ == '__main__':

    snp_l = []
    snp_cand_hc = read_candidates('hc', 'snp')
    snp_cand_fb = read_candidates('fb', 'snp')
    hc_snp_set = set(snp_cand_hc.v_id)
    fb_snp_set = set(snp_cand_fb.v_id)
    snp_cand = pandas.concat([snp_cand_hc, snp_cand_fb])
    snp_cand = snp_cand[~snp_cand.v_id.duplicated()]
    snp_cand['hc'] = snp_cand.v_id.apply(lambda i: i in hc_snp_set)
    snp_cand['fb'] = snp_cand.v_id.apply(lambda i: i in fb_snp_set)
    snp_cand['callers'] = snp_cand.apply(lambda row: concatClr(row), axis=1)

    indel_l = []
    indel_cand_hc = read_candidates('hc', 'indel')
    indel_cand_fb = read_candidates('fb', 'indel')
    hc_indel_set = set(indel_cand_hc.v_id)
    fb_indel_set = set(indel_cand_fb.v_id)
    indel_cand = pandas.concat([indel_cand_hc, indel_cand_fb])
    indel_cand = indel_cand[~indel_cand.v_id.duplicated()]
    indel_cand['hc'] = indel_cand.v_id.apply(lambda i: i in hc_indel_set)
    indel_cand['fb'] = indel_cand.v_id.apply(lambda i: i in fb_indel_set)
    indel_cand['callers'] = indel_cand.apply(lambda row: concatClr(row), axis=1)
