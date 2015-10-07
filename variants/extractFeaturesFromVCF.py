import variants
import ped
import func

variants = reload(variants)
ped = reload(ped)

infile_ped = '/mnt/scratch/asalomatov/data/SSCped/SSC.ped'
myped = ped.Ped(infile_ped, ['collection'])
myped.getParents(11006)
myped.addVcf()

infile_ped = '/mnt/scratch/asalomatov/data/columbia/pcgc_ped.txt'
myped = ped.Ped(infile_ped, [])
myped.getParents('1-00034')
myped.getFather('1-00034')
myped.
myped.ped.head()
myped.ped.shape
myped.addVcf(file_pat = '/mnt/scratch/asalomatov/data/columbia/vcf/deco/%s_%s-02_%s-01.annotated-deco.vcf.gz')
sum(myped.ped.vcf.notnull())

infile_vcf = '/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/rerun200fam/11006-HC-pm50-ann.vcf.gz'
#infile_vcf = '/mnt/scratch/asalomatov/data/columbia/vcf/deco/1-03173_1-03173-02_1-03173-01.annotated-deco.vcf.gz'
'/mnt/scratch/asalomatov/data/columbia/vcf/'
myvars = variants.Variants(infile_vcf, '11006', myped.ped)
myvars.readFromVcf()



myvars.samples
record = myvars.vcf_reader.next()
record.samples
myvars._colNamesFormat()

