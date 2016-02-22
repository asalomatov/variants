for i in test_sets:                                                                                                                                                          
    df = pandas.read_csv(i, usecols=[0,1], sep="\t")
    df['fam']=os.path.basename(i).split('.')[0]
    df[['fam','CHROM','POS']].to_csv('/mnt/scratch/asalomatov/data/SSC/wes/dnm_files/test/'+os.path.basename(i)+'.csv', header=False, index=False)
