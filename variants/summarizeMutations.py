#!/usr/bin/python

import pandas as pd
import sys, os

#indir, outdir = sys.argv[1:]

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




