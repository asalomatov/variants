import pandas
import sys, os
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import train
import numpy
import pandas
import ped


work_dirs = sys.argv[2:]
print work_dirs
output_file = sys.argv[1]

cum_metr_df = pandas.DataFrame()

pred_list = ['GBM_lvl9_stdFalse_cut0.5_splt0.99_300__0502_100_1_0.1_tstlvl10.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_500__0502_100_1_0.1_tstlvl10.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_1000__0502_100_1_0.1_tstlvl10.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_2000__0502_100_1_0.1_tstlvl10.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_5000__0502_100_1_0.1_tstlvl10.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_300__0502_tstlvl10.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_500__0502_tstlvl10.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_1000__0502_tstlvl10.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_2000__0502_tstlvl10.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_5000__0502_tstlvl10.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_300__0502_linear_1_balanced_tstlvl10.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_500__0502_linear_1_balanced_tstlvl10.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_1000__0502_linear_1_balanced_tstlvl10.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_2000__0502_linear_1_balanced_tstlvl10.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_5000__0502_linear_1_balanced_tstlvl10.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_300__0502_100_1_0.1_tstlvl9.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_500__0502_100_1_0.1_tstlvl9.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_1000__0502_100_1_0.1_tstlvl9.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_2000__0502_100_1_0.1_tstlvl9.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_5000__0502_100_1_0.1_tstlvl9.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_300__0502_tstlvl9.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_500__0502_tstlvl9.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_1000__0502_tstlvl9.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_2000__0502_tstlvl9.csv',
             'LogReg_lvl9_stdFalse_cut0.5_splt0.99_5000__0502_tstlvl9.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_300__0502_linear_1_balanced_tstlvl9.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_500__0502_linear_1_balanced_tstlvl9.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_1000__0502_linear_1_balanced_tstlvl9.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_2000__0502_linear_1_balanced_tstlvl9.csv',
             'SVM_lvl9_stdFalse_cut0.5_splt0.99_5000__0502_linear_1_balanced_tstlvl9.csv']             

all_pred_df = pandas.DataFrame()
for d in work_dirs:
    for p in pred_list:
        tst = train.TrainTest('x',
                              '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt',
                              ['status'],
                              ['descr'])
        my_pred = os.path.join(d, p)
        if not os.path.isfile(my_pred):
            print(my_pred)
            continue
        df = pandas.read_csv(my_pred)
        df = df[~df.test_var_id.duplicated()]
        all_pred_df = pandas.concat([all_pred_df, df])
        tst.test_set_y = df['test_labels']
        tst.pred_y = df['pred_labels']
        # tst.pred_y = (df['pred_prob'] > cut_off).astype(int)
        tst.pred_y_prob = df['pred_prob']
        tst.getMetrics()
        metr_df = tst.perf_mertics
        metr_df['method'] = my_pred
        cum_metr_df = pandas.concat([cum_metr_df, metr_df])

print cum_metr_df
cum_metr_df.to_csv(output_file, index=False)
