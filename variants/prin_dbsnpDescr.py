import variants



### dbSNP fields
myfile = '/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz'
myvars = variants.Variants(myfile, 'aaa')
myvars.describeInfoFields()
