from scipy import stats, integrate 
import matplotlib.pyplot as plt
import seaborn as sns

#plot univariate distributions
x = np.random.normal(size=1000)
y = np.random.normal(size=1000)
plt.hold(False)
plt.scatter(x, y)
plt.plot(x)
sns.distplot(x)
plt.clf()


# read a vcf file using PyVCF
who
import vcf

#'ios_mut_status_norm.vcf.gz'
infile = '/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/rerun200fam/11006-HC-pm50-ann.vcf.gz'
vcf_reader = vcf.Reader(open(infile, 'r'))
for i, record in enumerate(vcf_reader):
    print record.POS, record.INFO['DP']
#    if i > 2:
#        sys.exit()

vcf_reader = vcf.Reader(open(infile, 'r'))
vars_MT = vcf_reader.fetch()
vcf_reader == vars_MT
print vcf_reader.next()
print vcf_reader.next()
print vars_MT.next()
type(vars_MT)
type(vcf_reader)
record = vcf_reader.next()
type(record)
record.REF
record.ALT
record.ALT[0]
record.INFO
record.FORMAT['GT']
record.var_type
record.is_deletion
record.genotype(vcf_reader.samples[0]).gt_type
record.genotype(vcf_reader.samples[0]).gt_type
record.genotype(vcf_reader.samples[0]).site.ALT
vcf_reader.metadata
vcf_reader.infos
vcf_reader.contigs
#        if not append:
#            self.variants = pd.DataFrame()
#        if not chrom is None:
#            self.vcf_reader = self.vcf_reader.fetch(chrom, start, end)
#        if not self.samples is None and self.samples != self.vcf_reader.samples:
#            print >> sys.stderr, 'old file: ', 
#            raise Exception('Trying to append variants from different samples')
#

# read a vcf file using pybedtools
import pybedtools
a = pybedtools.BedTool('/mnt/ceph/asalomatov/SSC_Eichler/rerun/ssc11006/11006-JHC-vars.vcf.gz')
b = pybedtools.BedTool('/mnt/scratch/asalomatov/data/b37/b37.exome-pm50.bed')
a.head()
a.count()
b.head()
b.ann
a_and_b = a.intersect(b, u=True)
a_and_b.head()
a_and_b.count()
type(a_and_b)
a_and_b.c
#access intervals/lines
a.file_type
a
feature = a[0:3]
for f in feature:
    print f.chrom
feature = a[0:3]
f = feature.next()
f.chrom
f.fields
f.attrs

b[0]
aa = feature.next()
aa



# read Eichlers summarised denovos from Ioss and Crumm
import pandas as pd
import sys, os

fam_quad = ['13188', '14011', '11964', '13048', '11491', '13793', '11190', '13890', '13835', '12810', '12390', '13169', '12905', '11569', '11629', '11469', '12106', '11773', '13447', '12161', '13116', '11013', '11872', '11172', '11711', '11715', '12011', '14201', '12741', '11390', '11959', '13926', '13335', '11942', '13815', '12373', '12285', '13593', '12703', '11029', '11659', '11472', '11459', '11610', '11788', '13606', '11229', '13346', '11452', '11479', '11722', '13629', '12152', '12153', '12630', '12578', '11696', '12304', '13533', '12358', '12233', '11691']

fam_trio = ['11193', '11195', '11198', '11827', '13415', '11989', '13733', '11055', '11056', '11545', '11303', '12073', '12521', '11660', '11388', '11262', '11707', '13008', '12933', '13844', '11184', '11834', '12437', '12430', '11109', '12532', '11023', '11375', '13314', '13557', '13158', '12300', '11471', '13494', '13857', '12381', '11205', '13914', '13757', '12015', '13610', '14292', '12157', '13863', '13678', '11120', '13530', '13532', '11124', '12641', '11083', '11218', '13668', '13742', '11518', '13741', '13333', '12249', '11009', '11510', '12086', '12674', '11599', '13031', '11096', '11948', '11093', '11947', '11556', '11346', '11224', '13207', '12444', '11506', '11504', '12036', '11587', '12237', '12335', '12130', '11425', '12238', '14020', '12621', '13517', '11753', '12185', '11006', '11069', '11141', '12744', '11064', '11148', '11734', '11863', '12225', '12341', '12346', '12198', '11526', '11523', '13812', '11480', '11928', '12114', '12118', '11246', '12752', '12296', '12212', '14006', '11498', '11043', '12555', '12667', '13822', '12603', '11396', '11257', '13701', '11398', '13274', '11653', '11843', '11969']

families = fam_trio + fam_quad
print families
# read data into dataframes for these families

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

#first var in 11006 HC
#1       866511  rs60722469      C       CCCCT   1382.09 PASS   
