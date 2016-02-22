import ped
import sys

infile_ped = '/mnt/scratch/asalomatov/data/SSC/SSCped/SSC.ped'
myped = ped.Ped(infile_ped, ['collection'])
myped.addBam(file_pat='/mnt/scratch/asalomatov/data/SSC/wes/bam/%s.realigned.recal.bam')
myped.ped.dropna(subset=['bam'], inplace=True)
myped.ped.shape
myped.addBai(file_pat='/mnt/scratch/asalomatov/data/SSC/wes/bam/%s.realigned.recal.bam.bai')
myped.ped.dropna(subset=['bai'], inplace=True)
myped.ped.shape
myped.ped.head()

probands = myped.getAllProbands()
siblings = list(myped.ped.ind_id[myped.ped.ind_id.str.contains('s1')])
print 'probands ', len(probands)
print 'siblings ', len(siblings)

out_dir = '/mnt/scratch/asalomatov/data/SSC/wes/dnm_files/bam'
for p in probands + siblings:
    fa = myped.getChildsFather(p)
    mo = myped.getChildsMother(p)
    df = myped.ped[['ind_id', 'bam']][myped.ped.ind_id.isin([fa, mo, p])]
    print df.shape
    if df.shape[0] != 3:
        continue
    df.to_csv
    print df
    print ' '

