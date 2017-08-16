from __future__ import print_function
import pandas
import sys, os
import train

cut_off = int(sys.argv[1])
output_file = sys.argv[2]
pred_files = sys.argv[3:]
pred_files = [i for i in pred_files if 'vote' not in i]

mthd_thresholds = {'LogReg': 0.7,
                   'GBM': 0.3,
                   'SVM': 0.4}


def whichThre(model_name, thre_dict):
    for i in thre_dict.keys():
        if i in model_name:
            return thre_dict[i]
    return None

print(pred_files)
print(cut_off)
print(output_file)

all_pred_df = pandas.DataFrame()
for my_pred in pred_files:
    tst = train.TrainTest('x',
                          '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt',
                          ['status'],
                          ['descr'])
    if not os.path.isfile(my_pred):
        print(my_pred)
        continue
    df = pandas.read_csv(my_pred)
    df = df[~df.test_var_id.duplicated()]
    mthd_thr = whichThre(my_pred, mthd_thresholds)
    print('method cut off is %s' % mthd_thr)
    df.ix[:, 'pred_labels'] = (df['pred_prob'] > mthd_thr).astype(int)
    if 'GBM' in my_pred:
        df.ix[:, 'pred_labels'] = df.ix[:, 'pred_labels'] * 1
    print(df.pred_labels.value_counts())
    all_pred_df = pandas.concat([all_pred_df, df])

# aggregate by mutation id
x = all_pred_df.groupby('test_var_id').agg(
    {'DP_offspring': lambda x: x.iloc[0],
     'DP_father': lambda x: x.iloc[0],
     'DP_mother': lambda x: x.iloc[0],
     'pred_labels': 'sum',
     'test_labels': lambda x: x.iloc[0],
     'test_var_alleles': lambda x: x.iloc[0]}).reset_index()

tst.test_set_y = x['test_labels']
#tst.pred_y = x['pred_labels']
majority_pred = (x['pred_labels'] >= cut_off).astype(int)
tst.pred_y = majority_pred
#tst.pred_y_prob = x['pred_prob']
tst.getMetrics()
metr_x = tst.perf_mertics
metr_x['method'] = 'vote'

x[majority_pred == 1].to_csv(output_file,  index=False)
metr_x.to_csv(output_file+'.metrics', index=False)


