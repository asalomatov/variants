import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import variants
import func
import pandas
import os
import features

### SSC ped
#ped = reload(ped)
#func = reload(func)
infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
myped = ped.Ped(infile_ped, ['collection'])
myped.getParents(11006)
#myped.addVcf(file_pat = '/mnt/scratch/asalomatov/data/SSC/wes/vcf/raw/%s.family.vqsr.sorted.vcf.gz')
myped.addVcf()
myped.ped.dropna(subset=['vcf'], inplace=True)
myped.addBam()
myped.ped.dropna(subset=['bam'], inplace=True)


func = reload(func)
features = reload(features)
variants = reload(variants)

f = features.Features(myped, '/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_denovo_clean_snp.tsv')
f.initTrioFor('11006.p1')
f.father_bam
f.sample_vcf
f.extractFeatures()

myped.getIndivBAM('11006.p1')

# 

