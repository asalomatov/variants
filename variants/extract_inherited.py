#!/usr/bin/env python
from __future__ import print_function
import sys
import func
import features
import variants
import numpy
import pandas
import ped
import os
import yaml
import argparse
# from multiprocessing import Pool
import pkg_resources

arg_parser = argparse.ArgumentParser(
    description='Extract rare inherited mutations from VCF')
arg_parser.add_argument('child_id',
                        # nargs='+',
                        type=str,
                        help='A string identifying the child in a trio, must be\
                        present in the ped file, as well as in\
                        the corresponding vcf file headers')
arg_parser.add_argument('yaml_config_file',
                        # nargs='+',
                        type=str,
                        help='A config file, defining necessary variables.')
arg_parser.add_argument('--remove-tempfiles',
                        type=str,
                        default='yes',
                        help='yes/no')

args = arg_parser.parse_args()
print(args)
child_id = args.child_id
config_file = args.yaml_config_file
rm_tmp = (args.remove_tempfiles.lower() == 'yes')
print('remove_tempfiles is set to %s' % str(rm_tmp))

# get parameters from yaml config file
with open(config_file, 'r') as f:
    cfg = yaml.safe_load(f)
min_DP = cfg['min_DP']
var_type = cfg['variant_type']
ped_file_extended = cfg['ped_file_extended']
genome_ref = cfg['genome_ref']
output_dir = cfg['output_directory']
population_AF_max = cfg['population_AF_max']
population_AF_min = cfg['population_AF_min']

# create output dirs
func.makeDir(output_dir)
# populate ped DF
myped = ped.Ped(ped_file_extended, cfg['ped_extra_clmns'])

f = features.Features(myped, None)

# trio has to be complete with no files missing
if not f.initTrioFor(child_id):
    sys.stderr.write('\nfailed to initialize trio for ' + child_id)
    sys.exit(1)

sys.stdout.write('\ninitialized trio for ' + child_id)
sys.stdout.write('\n')
print('father and mother:')
print((f.father_id, f.mother_id))
# extract just the trio, drop non variant loci
# annotate, and filter based on AF
vrs = variants.Variants(f.sample_vcf, f.family_id)
vrs.readVcfToDF(sample_list=[f.sample_id, f.father_id, f.mother_id])
vrs.removeNoGT(f.sample_id)
vrs.removeHomRef(f.sample_id)
vrs.saveAsVcf(os.path.join(output_dir, 'aaa.vcf'))
