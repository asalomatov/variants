import pandas
import sys, os
import func

pred_file = sys.argv[1]
tag = sys.argv[2]

kv_vcf = pandas.read_csv('/mnt/scratch/asalomatov/data/columbia/feature_sets/known/all_known.txt', sep='\t')
kv_vcf = kv_vcf[['ind_id','CHROM', 'POS', 'REF_offspring', 'ALT_base_offspring', 'status', 'descr', 'DP_offspring', 'DP_father', 'DP_mother']]
kv_vcf = kv_vcf[kv_vcf.descr.isin(['after'])]
kv_vcf['var_id'] = kv_vcf.ind_id.astype(str)+'_'+kv_vcf.CHROM.astype(str)+'_'+kv_vcf.POS.astype(str)

mypred = pandas.read_csv(pred_file)
mypred['var_id'] = mypred['test_var_id']
mypred_u = mypred[~mypred.var_id.duplicated()]
mypred_u = mypred_u.merge(kv_vcf[['var_id', 'status']], on='var_id', how='left')
print 'status', mypred_u.status.value_counts()
print 'pred_labels', mypred_u.pred_labels.value_counts()
print 'test.labels', mypred_u.test_labels.value_counts()
c_status_known = ~mypred_u.status.isnull()
c_pred_pos = mypred_u.pred_labels == 1
c_status_pos = mypred_u.test_labels == 1
c_status_neg = mypred_u.test_labels == 0
mypred_u_res = mypred_u[c_status_known | c_pred_pos] 
print 'shape', mypred_u_res.shape
print 'status', mypred_u_res.status.value_counts()
print 'pred_labels', mypred_u_res.pred_labels.value_counts()
print 'test_labels', mypred_u_res.test_labels.value_counts()
func.writePredAsVcf(mypred_u_res,
                    '/mnt/scratch/asalomatov/data/columbia/vars_of_intr_gbm' + tag + '.tsv')

