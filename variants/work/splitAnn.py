# split per transcript annotations from tsv file generated by snpSift

import pandas
import sys

inp_file = sys.argv[1]

v = pandas.read_table(inp_file)
print v.shape
print v.columns
v.columns = [u'CHROM', u'POS', u'ID', u'REF', u'ALT', u'QUAL', u'VARTYPE',
       u'EFFECT', u'IMPACT', u'GENE', u'FEATURE',
       u'FEATUREID', u'BIOTYPE', u'ind_id', u'dbNSFP_genename',
       u'dbNSFP_Ensembl_transcriptid', u'dbNSFP_Ensembl_proteinid',
       u'dbNSFP_rs_dbSNP146', u'dbNSFP_aapos', u'dbNSFP_aaref',
       u'dbNSFP_aaalt', u'dbNSFP_Uniprot_acc_Polyphen2',
       u'dbNSFP_Uniprot_id_Polyphen2', u'dbNSFP_Uniprot_aapos_Polyphen2',
       u'dbNSFP_Ensembl_geneid', u'spidex_gene', u'spidex_strand',
       u'spidex_transcript', u'spidex_exon_number', u'spidex_location',
       u'spidex_cds_type', u'spidex_ss_dist']
print v.columns
v=v[~v.dbNSFP_aapos.isnull()]
print v.shape
v=v[v.dbNSFP_Ensembl_transcriptid!='.']
print v.shape
naa = v.dbNSFP_aapos.apply(lambda i: len(i.split(',')))
ntr = v.dbNSFP_Ensembl_transcriptid.apply(lambda i: len(i.split(',')))
v = v[ntr == naa]
print v.shape


def splitEff(row):
    eff = row['EFFECT'].split(',')
    print eff
    tr = row['FEATUREID'].split(',')
    print tr
    ind_id = row['ind_id']
    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    vartype = row['VARTYPE']
#    print aaref
    res = pandas.DataFrame({'ensemblTran': tr,
                            'effect': eff,
#                            'ensemblGene': ge,
                            'ind_id': [ind_id] * len(tr),
                            'chrom': [chrom] * len(tr),
                            'pos': [pos] * len(tr),
                            'ref': [ref] * len(tr),
                            'alt': [alt] * len(tr),
                            'vartype': [vartype] * len(tr)})
    return res
my_l = []
for i, row in v.iterrows():
    print i, row
    my_l.append(splitEff(row))

df_eff = pandas.concat(my_l)
df_miss = df_eff[df_eff.effect.str.contains('missense')]
df_miss['var_id'] = df_miss.ind_id + '_' + df_miss.chrom+'_'+ df_miss.pos.astype(str)+'_'+df_miss.ensemblTran

def splitCols(row):
    tr = row['dbNSFP_Ensembl_transcriptid'].split(',')
    print tr
    pr = row['dbNSFP_Ensembl_proteinid'].split(',')
    print pr
#    ge = row['dbNSFP_Ensembl_geneid'].split(',')
#    if len(ge) == 1:
#        ge = ge * len(pr)
#    print ge
    aapos = row['dbNSFP_aapos'].split(',')
    print aapos
    aaref = row['dbNSFP_aaref']
    print [aaref] * len(aapos)
    aaalt = row['dbNSFP_aaalt']
    ind_id = row['ind_id']
    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    vartype = row['VARTYPE']
#    print aaref
    res = pandas.DataFrame({'ensemblTran': tr,
                            'ensemblProt': pr,
#                            'ensemblGene': ge,
                            'aapos': aapos,
                            'aaref': [aaref] * len(aapos),
                            'aaalt': [aaalt] * len(aapos),
                            'ind_id': [ind_id] * len(aapos),
                            'chrom': [chrom] * len(aapos),
                            'pos': [pos] * len(aapos),
                            'ref': [ref] * len(aapos),
                            'alt': [alt] * len(aapos),
                            'vartype': [vartype] * len(aapos)})
    return res

res_l = []
for i, row in v.iterrows():
    print i, row
    res_l.append(splitCols(row))

df = pandas.concat(res_l)

df = df[['ind_id', 'chrom', 'pos', 'ref', 'alt', 'ensemblTran', 'ensemblProt', 'aapos', 'aaref', 'aaalt', 'vartype']]
df['var_id'] = df.ind_id+'_'+df.chrom+'_'+ df.pos.astype(str)+'_'+df.ensemblTran




#x = v.head().apply(splitCols, axis=1)
