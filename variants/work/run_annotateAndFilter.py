import os

clrs = ['fb', 'hc', 'pl']
btch = ['b1-2', 'b3', 'b4', 'b5']
btchyml = ['b_1-2', 'b_3', 'b_4', 'b_5']
vartype = ['INDEL', 'SNP']

cmd_l = []
for clr in clrs:
    for b in btch:
        byml = '_'.join([b[0], b[1:]])
        for vty in vartype:
            print clr, b, byml, vty
            if vty == 'INDEL':
                mdir = '/mnt/xfs1/home/asalomatov/projects/spark/denovo/indels'
            else:
                mdir = '/mnt/xfs1/home/asalomatov/projects/spark/denovo/snp'
            cmd1 = 'annotateAndFilter.py'
            cmd2 = os.path.join(mdir,
                                clr,
                                b,
                                clr + '-' + vty + '-class.csv')
            cmd3 = 'config_files/cfg_spark_'+byml+'_'+vty.lower()+'_'+clr+'.yml .3'
            cmd = ' '.join([cmd1, cmd2, cmd3])
            print cmd
            cmd_l.append(cmd)
