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

pred_list = ['GBM_lvl1_stdFalse_cut0.5_splt0.9_5000__100_1_0.1.csv',
             'GBM_lvl2_stdFalse_cut0.5_splt0.9_5000__100_1_0.1.csv',
             'GBM_lvl3_stdFalse_cut0.5_splt0.9_5000__100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_1000_ZZZ_500_4_0.2.pkl',
             'GBM_lvl4_stdFalse_cut0.5_splt0.9_5000__100_1_0.1.csv',
             'GBM_lvl6_stdFalse_cut0.5_splt0.9_5000__100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.9_5000__100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_1000_ZZZ_500_4_0.2.pklB']

#pred_list = ['SVM_lvl1_stdTrue_cut0.5_splt0.9_5000__linear_1_balanced.csv',
#             'SVM_lvl2_stdTrue_cut0.5_splt0.9_5000__linear_1_balanced.csv',
#             'SVM_lvl3_stdTrue_cut0.5_splt0.9_5000__linear_1_balanced.csv',
#             'SVM_lvl4_stdTrue_cut0.5_splt0.9_5000__linear_1_balanced.csv',
#             'SVM_lvl6_stdTrue_cut0.5_splt0.9_5000__linear_1_balanced.csv',
#             'SVM_lvl7_stdTrue_cut0.5_splt0.9_5000__linear_1_balanced.csv',
#             'SVM_lvl7_stdTrue_cut0.5_splt0.9_10000__linear_1_balanced.csv']

pred_list = ['GBM_lvl7_stdFalse_cut0.5_splt0.99_300_train01_0329_100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.99_500_train01_0329_100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.99_1000_train01_0329_100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.99_1500_train01_0329_100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.99_2000_train01_0329_100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.99_3000_train01_0329_100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.99_4000_train01_0329_100_1_0.1.csv',
             'GBM_lvl7_stdFalse_cut0.5_splt0.99_5000_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_300_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_500_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_1000_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_1500_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_2000_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_3000_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_4000_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_5000_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_300_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_500_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_1000_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_1500_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_2000_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_3000_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_4000_train01_0329_100_1_0.1.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_5000_train01_0329_100_1_0.1.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_300_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_500_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_1000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_1500_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_2000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_3000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_4000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl8_stdFalse_cut0.5_splt0.99_5000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_300_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_500_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_1000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_1500_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_2000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_3000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_4000_train01_0329_100_1_0.1_tstlvl7.csv',
             'GBM_lvl9_stdFalse_cut0.5_splt0.99_5000_train01_0329_100_1_0.1_tstlvl7.csv']

pred_list = ['dl_seq_lvl9_stdFalse_cut0.5_splt0.99_30000__0405_tstlvl10.csv',
             'dl_seq_lvl9_stdFalse_cut0.5_splt0.99_30000__0405_tstlvl9.csv']
             
for d in work_dirs:
    for p in pred_list:
        tst = train.TrainTest('x',
                              '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt',
                              ['status'],
                              ['descr'])
        my_pred = os.path.join(d, p)
        if not os.path.isfile(my_pred):
            continue
        df = pandas.read_csv(my_pred)
        df = df[~df.test_var_id.duplicated()]
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
