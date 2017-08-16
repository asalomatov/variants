'''
Output a table with mutaions ready for validation submittion.
'''

from collections import OrderedDict
import pandas
import func


validation_clms = OrderedDict()
validation_clms['SP_id'] = 'Proband SPID'
validation_clms['SP_id_mo'] = 'Mother SPID'
validation_clms['SP_id_fa'] = 'Father SPID'
validation_clms['CHROM'] = 'Chromosome'
validation_clms['POS'] = 'Chromosome position'
# validation_clms['POS_end'] = 'Chromosome End'
validation_clms['REF'] = 'Reference Allele'
validation_clms['ALT'] = 'Alternative Allele'
validation_clms['Zygosity'] = 'Zygosity'
validation_clms['VARTYPE'] = 'Type of mutation'
validation_clms['ANN.GENE'] = 'Gene name'
validation_clms['Feature'] = 'Transcript ID'
validation_clms['ANN.HGVS_C'] = 'HGVSc'
validation_clms['ENSP'] = 'Protein ID'
# validation_clms['ANN.RANK'] = 'Exon'
validation_clms['ANN.HGVS_P'] = 'HGVSp'
validation_clms['ANN.HGVS_P'] = 'HGVSp'
# validation_clms['SP_id_sibs'] = 'Siblings'
validation_clms['Siblings'] = 'Siblings'
validation_clms['Consented siblings'] = 'Consented siblings'
validation_clms['Siblings with dna'] = 'Siblings with dna'
validation_clms['Sequenced siblings'] = 'Sequenced siblings'
validation_clms['Affected siblings'] = 'Affected siblings'
validation_clms['Unaffected siblings'] = 'Unaffected siblings'
validation_clms['Affected twins'] = 'Affected twins'
validation_clms['Unaffected twins'] = 'Unaffected twins'
validation_clms['HGVSc;Exon;Intron;HGVSp'] = 'HGVSc;Exon;Intron;HGVSp'

ped_file = '/mnt/xfs1/scratch/asalomatov/data/SPARK/ped/b1-9_pedigree.ped'
dmrt_df = pandas.read_csv('/mnt/xfs1/scratch/asalomatov/data/SPARK/info/SPARK_datamart_validation.csv')
AAA = pandas.read_csv(
    '/mnt/ceph/users/asalomatov/spark/denovo/def1/spark_b1-9_for_validation-1_LOF_DMIS_ALL_clinical.csv')
if 'SP_id' not in AAA.columns:
    AAA['SP_id'] = AAA.ind_id
AAA['Zygosity'] = 'Het'
AAA['SP_id_fa'] = AAA['SP_id'].apply(
    lambda i: func.parentForOffspring(i, ped_file, sex=1))
AAA['SP_id_mo'] = AAA['SP_id'].apply(
    lambda i: func.parentForOffspring(i, ped_file, sex=2))
# AAA['SP_id_sibs'] = AAA['SP_id'].apply(
#     lambda i: func.siblingsForOffspring(i, ped_file))
dmrt_columns = AAA.SP_id.apply(lambda i: func.sibFromDataMart(i, dmrt_df))
AAA = AAA.merge(dmrt_columns, how='left', left_index=True, right_index=True)

AAA = AAA[validation_clms.keys()]
AAA.columns = [validation_clms[i] for i in AAA.columns]
AAA['Lab Name'] = None
AAA['Lab ID'] = None
AAA['Date Collected'] = None
AAA['Date Received'] = None
AAA['Date Reported'] = None
AAA['Date Ordered'] = None
