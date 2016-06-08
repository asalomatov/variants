from __future__ import print_function
import sys
import func
import train
import numpy
import pandas
import ped
import os
from sklearn.externals import joblib
#from keras.models import model_from_json
import yaml
import argparse

arg_parser = argparse.ArgumentParser(
    description='Classify De Novo mutations.')
arg_parser.add_argument('child_id',
                        nargs='+',
                        type=str,
                        help='A string identifying the child in a trio, must be\
                        present in the ped file, as well as in the corresponding\
                        vcf file headers')
arg_parser.add_argument('yaml_config_file',
                        nargs='+',
                        type=str,
                        help='A config file, defining necessary variables.')
arg_parser.add_argument('class_probability_threshold',
                        nargs='+',
                        type=float,
                        help='Only mutations with scores higher than this\
                        will be found in the output.')
arg_parser.add_argument('--sklearn_model_pkl',
                        nargs='+',
                        type=str,
                        default=' ',
                        help='A pickle of a classifier. If not a file, an internal\
                        model will be used.')

args = arg_parser.parse_args()
print(args)
child_id = args.child_id[0]
config_file = args.yaml_config_file[0]
prob_cutoff = args.class_probability_threshold[0]
print(prob_cutoff)

model = args.sklearn_model_pkl[0]
print(model)

sys.exit(1)

myped = ped.Ped(args.ped_file[0])
myped.addBam(field='ind_id', file_pat=args.bam_pattern[0])



child_id = sys.argv[1]
model = sys.argv[3]
lvl = int(sys.argv[4])
is_keras = bool(int(sys.argv[5]))
prob_cutoff = float(sys.argv[6])

# get parameters from yaml config file
with open(config_file, 'r') as f:
    cfg = yaml.safe_load(f)
min_DP = cfg['min_DP']
var_type = cfg['variant_type']
vcf_pat = cfg['vcf_pattern']
ped_file = cfg['ped_file']
known_vars = cfg['known_variants']
output_dir = cfg['output_directory']
test_set_pat = output_dir + '/%s'
m_pkl = joblib.load(model)
list_of_features = m_pkl['features']
if lvl == 0:
    m_pkl['extra_col_names'] = []
    m_pkl['y_name'] = []

myped = ped.Ped(ped_file)
myped.addTestFile(field='ind_id', file_pat=test_set_pat)
myped.ped.dropna(subset=['test'], inplace=True)
myped.ped.reset_index(inplace=True)

test_labels = numpy.array([], dtype=int)
pred_labels = numpy.array([], dtype=int)
test_var_id = numpy.array([], dtype=str)
test_alleles = numpy.array([], dtype=str)
pred_prob = numpy.array([], dtype=float)
dp_offspring = numpy.array([], dtype=int)
dp_father = numpy.array([], dtype=int)
dp_mother = numpy.array([], dtype=int)

for i, row in myped.ped.iterrows():
    if row['ind_id'] != child_id:
        continue
#    print 'processing', i, row['ind_id']
    #myped.ped.test.iat[i]
    tst = train.TrainTest(row['test'],
                          list_of_features,
                          m_pkl['y_name'],
                          m_pkl['extra_col_names'] + ['DP_offspring', 'DP_father', 'DP_mother'])
    if is_keras:
        tst.is_keras = True
    tst.feature_list = list_of_features
    tst.readDataSet()
#    print 'data_set shape is ', tst.data_set.shape
    if tst.data_set.empty:
        continue
    tst.addLabels(level=lvl)
#    print 'data_set shape is ', tst.data_set.shape
#    print tst.data_set.label.value_counts()
    tst.dropNA('label')
#    print 'data_set shape with non null labels is ', tst.data_set.shape
    if tst.is_keras:
        tst.model = model_from_json(m_pkl['model'])
        tst.model.load_weights(m_pkl['weights_file'])
    else:
        tst.model = m_pkl['model']
    tst.stdize = m_pkl['stdize']  # bool(int(sys.argv[6]))
    tst.trshold = m_pkl['threshold']  # float(sys.argv[7])
    tst.train_set_var_id = m_pkl['train_var_id']
    tst.data2Test()
#    print 'test_set_X shape is', tst.test_set_X.shape
    tst.predictClass(tst.threshold)
    #tst.getMetrics()
    test_labels = numpy.concatenate((test_labels, tst.test_set_y))
    pred_labels = numpy.concatenate((pred_labels, tst.pred_y))
    pred_prob = numpy.concatenate((pred_prob, tst.pred_y_prob))
    test_var_id = numpy.concatenate((test_var_id, tst.test_set_var_id))
    test_alleles = numpy.concatenate((test_alleles, tst.test_set_alleles))
    dp_offspring = numpy.concatenate((dp_offspring, tst.test_set_DP_offspring))
    dp_father = numpy.concatenate((dp_father, tst.test_set_DP_father))
    dp_mother = numpy.concatenate((dp_mother, tst.test_set_DP_mother))
res = pandas.DataFrame({'test_labels': test_labels, 'pred_labels': pred_labels,
                        'pred_prob': pred_prob, 'test_var_id': test_var_id,
                        'test_var_alleles': test_alleles,
                        'DP_offspring': dp_offspring,
                        'DP_father': dp_father,
                        'DP_mother': dp_mother})
m_name = os.path.basename(model)
m_name = '.'.join(m_name.split('.')[:-1]) + '_tstlvl' + str(lvl)
res['method'] = m_name
res = res[~res.test_var_alleles.str.contains('nan')]
res.to_csv(os.path.join(output_dir, m_name + '.csv'), index=False)
res['var_id'] = res['test_var_id']
res_u = res[~res.var_id.duplicated()]
res_u.reset_index(inplace=True)
res_u.ix[:, 'pred_labels'] = (res_u['pred_prob'] > prob_cutoff).astype(int)
#res_u = res_u[res_u.pred_labels == 1]
#outp_tsv = os.path.join(output_dir, m_name + '.tsv')
outp_tsv = os.path.join(output_dir, child_id + '.tsv')
func.writePredAsVcf(res_u, outp_tsv, min_DP=min_DP)

script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
cmd = ' '.join([os.path.join(script_dir, 'vcf2table.sh'),
               outp_tsv,
               script_dir,
               child_id])
print(cmd)
func.runInShell(cmd)

summarizeVariants.summarizeMutations(os.path.join(output_dir, child_id + '-ann-onePline.tsv'),
                                                  os.path.join(output_dir, 'denovo'),
                                                  config_file)


#cmd = ' '.join([os.path.join(script_dir, 'work', 'summarizeMutations.py'),
#               os.path.join(output_dir, child_id + '-ann-onePline.tsv'),
#               os.path.join(output_dir, 'denovo'),
#                config_file])
#func.runInShell(cmd)

# work/summarizeMutations.py /mnt/xfs1/home/asalomatov/projects/spark/feature_sets/hc/trio003.p1_642940-ann-onePline.tsv /mnt/xfs1/home/asalomatov/projects/spark/feature_sets/hc/denovo cfg_spark.yml
