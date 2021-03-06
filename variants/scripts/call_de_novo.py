#!/usr/bin/env python
from __future__ import print_function
import sys
from variants import func
from variants import features
from variants import train
import numpy
import pandas
from variants import ped
import os
from sklearn.externals import joblib
from keras.models import model_from_json
import yaml
import argparse
from multiprocessing import Pool
import pkg_resources


# wrap a funtion for use with multiprocessing
def multi_wrap_readBamReadcount(args):
    return func.readBamReadcount(*args)

arg_parser = argparse.ArgumentParser(
    description='Classify De Novo mutations.')
arg_parser.add_argument('child_id',
#                        nargs='+',
                        type=str,
                        help='A string identifying the child in a trio, must be\
                        present in the ped file, as well as in the corresponding\
                        vcf file headers')
arg_parser.add_argument('yaml_config_file',
#                        nargs='+',
                        type=str,
                        help='A config file, defining necessary variables.')
arg_parser.add_argument('class_probability_threshold',
#                        nargs='+',
                        type=float,
                        help='Only mutations with scores higher than this\
                        will be found in the output.')
arg_parser.add_argument('--sklearn_model_pkl',
#                        nargs='+',
                        type=str,
                        default=' ',
                        help='A pickle of a classifier. If not a file, an internal\
                        model will be used.')

arg_parser.add_argument('--remove_tempfiles',
                        type=str,
                        default='yes',
                        help='yes/no')

args = arg_parser.parse_args()
print(args)
child_id = args.child_id
config_file = args.yaml_config_file
prob_cutoff = args.class_probability_threshold
model = args.sklearn_model_pkl
rm_tmp = (args.remove_tempfiles.lower() == 'yes')
print('remove_tempfiles is set to %s' % str(rm_tmp))

# get parameters from yaml config file
with open(config_file, 'r') as f:
    cfg = yaml.safe_load(f)
min_DP = cfg['min_DP']
var_type = cfg['variant_type']
# vcf_pat = cfg['vcf_pattern']
# bam_pat = cfg['bam_pattern']
# bai_pat = bam_pat + '.bai'
# ped_file = cfg['ped_file']
ped_file_extended = cfg['ped_file_extended']
bam_readcount = cfg['bam_readcount']
genome_ref = cfg['genome_ref']

known_vars = None
output_dir = cfg['output_directory']
test_set_pat = output_dir + '/%s'
dnv_def = cfg['denovo_definition']

if not os.path.isfile(model):
    # script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
    # print(script_dir)
    model_dir = pkg_resources.resource_filename('variants',
                                                'denovo_classifier_model_' +
                                                var_type.upper())
    mdl_files = os.listdir(model_dir)
    model = [i for i in mdl_files if '.pkl_' not in i]
    if len(model) != 1:
        print('found %s models' % len(model))
        print('supply one model via --sklearn_model_pkl argument')
        sys.exit(1)
    model = os.path.join(model_dir, model[0])
    print(model)

m_pkl = joblib.load(model)
list_of_features = m_pkl['features']
# hardcode lvl, this is intended for external use only
lvl = 0
is_keras = bool(int(m_pkl['is_keras']))

if lvl == 0:
    m_pkl['extra_col_names'] = []
    m_pkl['y_name'] = []

# create output dirs
func.runInShell('mkdir -p ' + output_dir)
if known_vars:
    output_dir_known = os.path.join(output_dir, 'known')
    func.runInShell('mkdir -p ' + output_dir_known)

# populate ped DF
myped = ped.Ped(ped_file_extended, ['bam', 'vcf'])
f = features.Features(myped, known_vars)

# trio has to be complete with no files missing
if not f.initTrioFor(child_id):
    sys.stderr.write('\nfailed to initialize trio for ' + child_id)
    sys.exit(1)
else:
    sys.stdout.write('\ninitialized trio for ' + child_id)
    sys.stdout.write('\n')
    f.extractFeatures(genome_ref, bam_readcount, var_type, dnv_def)
    pool = Pool(3)
    fam_features = pool.map(multi_wrap_readBamReadcount,
                            [(f.sample_features, var_type, 14),
                             (f.father_features, var_type, 14),
                             (f.mother_features, var_type, 14)])
    if rm_tmp:
        f.removeTmpDir()
    fam_features[0].columns = ['CHROM', 'POS'] + func.addSuffix(
        fam_features[0].columns[2:], '_offspring')
    fam_features[1].columns = ['CHROM', 'POS'] + func.addSuffix(
        fam_features[1].columns[2:], '_father')
    fam_features[2].columns = ['CHROM', 'POS'] + func.addSuffix(
        fam_features[2].columns[2:], '_mother')
    fam_f = fam_features[0]
    fam_f = fam_f.merge(fam_features[1], how='inner', on=['CHROM', 'POS'])
    fam_f = fam_f.merge(fam_features[2], how='inner', on=['CHROM', 'POS'])
    fam_f = fam_f[fam_f.ALT_count_offspring.astype(int) != 0]  # remove not-variant loci
    fam_f['ind_id'] = child_id
    # mark verified if information is available
    if known_vars:
        f.verified_variants.columns = ['ind_id', 'CHROM', 'POS',
                                       'ref', 'alt','status',  'descr']
        # deletions should have POS + 1
        c1 = f.verified_variants.ref.apply(len) > 1  # deletion 1
        c2 = f.verified_variants.alt.apply(len) == 1  # deletion 2
        f.verified_variants.ix[:, 'POS'] = f.verified_variants.POS.astype(int)
        f.verified_variants.ix[c1 & c2, 'POS'] = f.verified_variants.POS[c1 & c2] + 1
        f.verified_variants['POS'] = f.verified_variants['POS'].astype(str)
        fam_f = fam_f.merge(f.verified_variants[['ind_id', 'CHROM', 'POS', 'status', 'descr']],
                            how='left', on=['ind_id', 'CHROM', 'POS'])
        fam_f.ix[fam_f.status.isnull(), 'status'] = 'NA'
    clmns_to_drop = ['REF_avg_se_mapping_quality_offspring',
                     'REF_avg_se_mapping_quality_father',
                     'REF_avg_se_mapping_quality_mother',
                     'ALT_avg_se_mapping_quality_offspring',
                     'ALT_avg_se_mapping_quality_father',
                     'ALT_avg_se_mapping_quality_mother',
                     'INDEL_avg_se_mapping_quality_offspring',
                     'INDEL_avg_se_mapping_quality_father',
                     'INDEL_avg_se_mapping_quality_mother',
                     # 'ALT_base_offspring',
                     # 'ALT_base_father',
                     # 'ALT_base_mother',
                     'INDEL_base_offspring',
                     'INDEL_base_father',
                     'INDEL_base_mother']
#    fam_f.drop(clmns_to_drop, axis=1, inplace=True)
    fam_f['DP_offspring'] = fam_f['REF_count_offspring'] + fam_f['ALT_count_offspring']
    fam_f['DP_father'] = fam_f['REF_count_father'] + fam_f['ALT_count_father']
    fam_f['DP_mother'] = fam_f['REF_count_mother'] + fam_f['ALT_count_mother']
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring)
    sum(fam_f.DP_offspring - fam_f.REF_count_offspring - fam_f.ALT_count_offspring != 0)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father)
    sum(fam_f.DP_father - fam_f.REF_count_father - fam_f.ALT_count_father != 0)
    # filter out DP < 8
    c1 = fam_f['DP_offspring'] > min_DP - 1
    c2 = fam_f['DP_father'] > min_DP - 1
    c3 = fam_f['DP_mother'] > min_DP - 1
    fam_f = fam_f[c1 & c2 & c3]
    fam_f.to_csv(os.path.join(output_dir, child_id), sep='\t', index=False)
    myped.addTestFile(field='ind_id', file_pat=test_set_pat)
    myped.ped.dropna(subset=['test'], inplace=True)
    myped.ped.reset_index(inplace=True)

    test_labels = numpy.array([], dtype=int)
    pred_labels = numpy.array([], dtype=int)
    test_var_id = numpy.array([], dtype=str)
    test_alleles = numpy.array([], dtype=str)
    test_alleles_fa = numpy.array([], dtype=str)
    test_alleles_mo = numpy.array([], dtype=str)
    pred_prob = numpy.array([], dtype=float)
    dp_of = numpy.array([], dtype=int)
    dp_fa = numpy.array([], dtype=int)
    dp_mo = numpy.array([], dtype=int)
    num_alt_of = numpy.array([], dtype=int)
    num_alt_fa = numpy.array([], dtype=int)
    num_alt_mo = numpy.array([], dtype=int)
    num_alt_all = numpy.array([], dtype=int)
    inherit_fa = numpy.array([], dtype=int)
    inherit_mo = numpy.array([], dtype=int)
    inherit_prnts = numpy.array([], dtype=int)
    for i, row in myped.ped.iterrows():
        if row['ind_id'] != child_id:
            continue
    #    print 'processing', i, row['ind_id']
        # myped.ped.test.iat[i]
        tst = train.TrainTest(row['test'],
                              list_of_features,
                              m_pkl['y_name'],
                              m_pkl['extra_col_names'])
# +\
#                              ['DP_offspring', 'DP_father', 'DP_mother'])
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
        tst.stdize = m_pkl['stdize']
        tst.trshold = prob_cutoff
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
        test_alleles_fa = numpy.concatenate((test_alleles_fa,
                                             tst.test_set_alleles_fa))
        test_alleles_mo = numpy.concatenate((test_alleles_mo,
                                             tst.test_set_alleles_mo))
        dp_of = numpy.concatenate((dp_of, tst.test_set_DP_of))
        dp_fa = numpy.concatenate((dp_fa, tst.test_set_DP_fa))
        dp_mo = numpy.concatenate((dp_mo, tst.test_set_DP_mo))
        num_alt_of = numpy.concatenate((num_alt_of, tst.test_set_num_alt))
        num_alt_fa = numpy.concatenate((num_alt_fa, tst.test_set_num_alt_fa))
        num_alt_mo = numpy.concatenate((num_alt_mo, tst.test_set_num_alt_mo))
        num_alt_all = numpy.concatenate((num_alt_all,
                                         tst.test_set_num_alt_all))
        inherit_fa = numpy.concatenate((inherit_fa, tst.test_set_inherit_fa))
        inherit_mo = numpy.concatenate((inherit_mo, tst.test_set_inherit_mo))
        inherit_prnts = numpy.concatenate((inherit_prnts,
                                           tst.test_set_inherit_prnts))
    res = pandas.DataFrame({'test_labels': test_labels,
                            'pred_labels': pred_labels,
                            'pred_prob': pred_prob,
                            'var_id': test_var_id,
                            'alleles_of': test_alleles,
                            'alleles_fa': test_alleles_fa,
                            'alleles_mo': test_alleles_mo,
                            'num_alt_of': num_alt_of,
                            'num_alt_fa': num_alt_fa,
                            'num_alt_mo': num_alt_mo,
                            'num_alt_all': num_alt_all,
                            'inherit_fa': inherit_fa,
                            'inherit_mo': inherit_mo,
                            'inherit_prnts': inherit_prnts,
                            'DP_of': dp_of,
                            'DP_fa': dp_fa,
                            'DP_mo': dp_mo})
    m_name = os.path.basename(model)
    m_name = '.'.join(m_name.split('.')[:-1]) + '_tstlvl' + str(lvl)
    res['method'] = m_name
    res = res[~res.alleles_of.str.contains('nan')]
    # res['var_id'] = res['test_var_id']
    res_u = res[~res.var_id.duplicated()]
    res_u.reset_index(inplace=True)
    res_u.ix[:, 'pred_labels'] = (res_u['pred_prob'] > prob_cutoff).astype(int)
    res_u = res_u[res_u.pred_labels == 1]
    res_u.reset_index(inplace=True)
    if res_u.empty:
        print('found no de novo mutations')
    else:
        print('found de novo')
        print(res_u.shape)
        varid = res_u.var_id.apply(func.splitVarId)
        lls = res_u.alleles_of.apply(func.splitAlleles)
        res_u = res_u.merge(varid, left_index=True, right_index=True)
        res_u = res_u.merge(lls, left_index=True, right_index=True)
        res_u['POS'] = res_u.POS.astype(int)
        if var_type.lower() == 'indel':
            c_ins = res_u.ALT.str.contains('+', regex=False)
            c_del = res_u.ALT.str.contains('-', regex=False)
            res_u.ix[:, 'ALT'] = res_u.ALT.apply(lambda x: x.strip('+'))
            res_u.ix[:, 'ALT'] = res_u.ALT.apply(lambda x: x.strip('-'))
            res_u.ix[c_ins, 'ALT'] = res_u.REF[c_ins] + res_u.ALT[c_ins]
            res_u.ix[c_del, 'REF'] = res_u.ALT[c_del]
            res_u.ix[c_del, 'POS'] -= 1
            ref_pos = res_u[c_del].apply(lambda row:
                                         func.refAtPos(row['CHROM'],
                                                       row['POS'],
                                                       genref=genome_ref),
                                         axis=1)
            res_u.ix[c_del, 'REF'] = ref_pos + res_u.ix[c_del, 'REF']
            res_u.ix[c_del, 'ALT'] = ref_pos
        res_u['var_id'] = res_u.ind_id.astype(str) + '_' +\
                          res_u.CHROM.astype(str) + '_' +\
                          res_u.POS.astype(str)
        res_u['var_id_a'] = res_u.ind_id.astype(str) + '_' +\
                            res_u.CHROM.astype(str) + '_' +\
                            res_u.POS.astype(str) + '_' +\
                            res_u.REF + '_' +\
                            res_u.ALT
        res_u.DP_of = res_u.DP_of.astype(float).astype(int)
        res_u.DP_fa = res_u.DP_fa.astype(float).astype(int)
        res_u.DP_mo = res_u.DP_mo.astype(float).astype(int)
        res_u.num_alt_all = res_u.num_alt_all.astype(float).astype(int)
        res_u[['ind_id',
               'CHROM',
               'POS',
               'REF',
               'ALT',
               'pred_prob',
               'DP_of',
               'DP_fa',
               'DP_mo',
               'var_id',
               'var_id_a',
               'alleles_of',
               'alleles_fa',
               'alleles_mo',
               'num_alt_of',
               'num_alt_fa',
               'num_alt_mo',
               'num_alt_all',
               'inherit_fa',
               'inherit_mo',
               'inherit_prnts']].to_csv(os.path.join(output_dir,
                                                     child_id +
                                                     '-' + var_type +
                                                     '-class.csv'),
                                        index=False)
