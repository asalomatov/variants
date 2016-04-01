import pandas
import sys, os
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import train
import numpy
import pandas
import ped
import features

output_file = sys.argv[1]
dnm_pred_file = sys.argv[2]
lvl = int(sys.argv[3])
print dnm_pred_file

known_vars = '/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/ssc_denovo_clean_snp.tsv'
cum_metr_df = pandas.DataFrame()

kv = pandas.read_csv(known_vars, sep='\t', dtype=str)
kv.columns = ['ind_id', 'CHROM', 'POS', 'ref', 'alt', 'status',  'descr', 'vartype']
print kv.head()

#f = features.Features(None, known_vars)
#f.verified_variants.columns = ['ind_id', 'CHROM', 'POS', 'ref', 'alt', 'status',  'descr']
#f.verified_variants['POS'] = f.verified_variants['POS'].astype(str)

cum_metr_df = pandas.DataFrame()

dnm = pandas.read_csv(dnm_pred_file, header=None, dtype=str)
dnm.columns = ['ind_id', 'fam_id', 'CHROM', 'POS', 'SCORE']
dnm = dnm.merge(kv[['ind_id', 'CHROM', 'POS', 'status', 'descr']],
                how='left', on=['ind_id', 'CHROM', 'POS'])
print dnm.head()
tst = train.TrainTest('x',
                      '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt',
                      ['status'],
                      ['descr'])
tst.data_set = dnm
tst.addLabels(level=lvl)
tst.dropNA('label')
tst.pred_y_prob = tst.data_set['SCORE'].astype(float)
tst.test_set_y = tst.data_set['label'].astype(int)

for cut_off in [0.5, 0.75, 0.9]:
    tst.pred_y = (tst.pred_y_prob > cut_off).astype(int)
    tst.getMetrics()
    metr_df = tst.perf_mertics
    metr_df['method'] = 'DNMfilter'
    metr_df['min_score'] = cut_off
    cum_metr_df = pandas.concat([cum_metr_df, metr_df])


print cum_metr_df
cum_metr_df.to_csv(output_file, index=False)
