#!/mnt/xfs1/home/asalomatov/miniconda2/bin/python
from __future__ import print_function
import sys
from variants import func
import pandas
import os
import yaml
from variants import summarizeVariants
import argparse
import pkg_resources
import tempfile


print('command line:')
print(sys.argv)

arg_parser = argparse.ArgumentParser(
    description='Annotate and filter de novo mutations')
arg_parser.add_argument('input_file',
                        nargs='+',
                        type=str,
                        help='A file that looks like output\
                        of call_de_novo.py')
arg_parser.add_argument('yaml_config_file',
                        nargs='+',
                        type=str,
                        help='A config file, defining necessary variables.')
arg_parser.add_argument('class_probability_threshold',
                        nargs='+',
                        type=float,
                        help='Only mutations with scores higher than this\
                        will be found in the output.')

args = arg_parser.parse_args()
print(args)
input_file = args.input_file[0]
config_file = args.yaml_config_file[0]
prob_cutoff = args.class_probability_threshold[0]



# get parameters from yaml config file
with open(config_file, 'r') as f:
    cfg = yaml.safe_load(f)
min_DP = cfg['min_DP']
var_type = cfg['variant_type']
vcf_pat = cfg['vcf_pattern']
ped_file = cfg['ped_file']
known_vars = cfg['known_variants']
output_dir = cfg['output_directory']

dnvo = pandas.read_csv(input_file)
dnvo.ix[:, 'pred_labels'] = (dnvo['pred_prob'] > prob_cutoff).astype(int)
dnvo = dnvo[dnvo.pred_labels == 1]
dnvo.reset_index(inplace=True)
if dnvo.empty:
    sys.exit('No mutations at score %s' % prob_cutoff)

tmp_dir = tempfile.mkdtemp()
print(tmp_dir)
input_file_bn = os.path.splitext(os.path.basename(input_file))[0]
outp_tsv = os.path.join(tmp_dir, input_file_bn + '.tsv')
print(outp_tsv)
func.writePredAsVcf(dnvo, outp_tsv, min_DP=min_DP)
# script_name = os.path.basename(os.path.realpath(sys.argv[0]))
script_name = os.path.abspath(pkg_resources.resource_filename('variants',
                                                              'vcf2table.sh'))
cmd = ' '.join([script_name,
               outp_tsv,
               os.path.dirname(script_name),
                input_file_bn])
print(cmd)
func.runInShell(cmd)
summarizeVariants.summarizeMutations(os.path.join(tmp_dir, input_file_bn +
                                                  '-ann-onePline.tsv'),
                                     input_file_bn,
                                     output_dir,
                                     config_file)
