data = pyensembl.EnsemblRelease(75)
data.reference_name
x = data.genes()
y = filter(lambda i: i.biotype in 'protein_coding', x)
len(y)
mygenes=[]
for g in y:
    mygenes.append([g.contig, str(g.start - 1), str(g.start - 1 + g.length), g.name, g.id])

len(mygenes)
mygenes[0]
'\t'.join(mygenes[0])
import sys
myfile = open('/mnt/scratch/asalomatov/data/b37/ensembl-genes-raw.bed', 'a')
myfile
for line in mygenes:
        myfile.write('\t'.join(line)+'\n')
myfile.close()
                                
