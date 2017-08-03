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
import tempfile
# from multiprocessing import Pool
import pkg_resources
import logging

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
arg_parser.add_argument('--annotate',
                        type=str,
                        default='yes',
                        help='yes/no. If no, vcf or vep file is assumed to be\
                        annotated witn AF. If yes, the file will be annotated\
                        first.')
arg_parser.add_argument('--pipeline_directory',
                        type=str,
                        default='/mnt/xfs1/home/asalomatov/projects/pipeline/ppln',
                        help='Path to the folder with pipeline')


logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
args = arg_parser.parse_args()
print(args)
child_id = args.child_id
config_file = args.yaml_config_file
rm_tmp = (args.remove_tempfiles.lower() == 'yes')
logging.info('remove_tempfiles is set to %s' % str(rm_tmp))

# get parameters from yaml config file
with open(config_file, 'r') as f:
    cfg = yaml.safe_load(f)
min_DP = cfg['min_DP']
var_type = cfg['variant_type']
ped_file_extended = cfg['ped_file_extended']
genome_ref = cfg['genome_ref']
targ_bed = cfg['target_bed']
output_dir = cfg['output_directory']
population_AF_max = cfg['population_AF_max']
population_AF_min = cfg['population_AF_min']
genome_build = int(cfg['genome_build'])
ppln_dir = args.pipeline_directory

if genome_build == 19 or genome_build == 37:
    incl_make = os.path.join(ppln_dir, 'include.mk')
elif genome_build == 38:
    incl_make = os.path.join(ppln_dir, 'include_hg38.mk')
else:
    logging.critical('invalid genome build specified, only\
    builds 19, 37, 38 are supported.')
    sys.exit('Only builds 19, 37, 38 are supported')

logging.info('Using %s as pipeline config file' % incl_make)

# create output dirs
func.makeDir(output_dir)
# create temp dir 
tmp_dir = tempfile.mkdtemp()
logging.info('working dir is %s' % tmp_dir)
# get path to script
script_name = os.path.abspath(
    pkg_resources.resource_filename('variants', 'vep_snpeff_for_vcf.sh'))

# populate ped DF
myped = ped.Ped(ped_file_extended, cfg['ped_extra_clmns'])
f = features.Features(myped, None)

# trio has to be complete with no file missing
if not f.initTrioFor(child_id):
    logging.critical('\nfailed to initialize trio for ' + child_id)
    sys.exit(1)

sys.stdout.write('\ninitialized trio for ' + child_id)
sys.stdout.write('\n')
logging.info('father and mother: ' +
             ' '.join([f.father_id, f.mother_id]))
# extract just the trio, drop non variant loci
# annotate, and filter based on AF
vrs = variants.Variants(f.sample_vcf, f.family_id)
vrs.readVcfToDF(sample_list=[f.sample_id, f.father_id, f.mother_id])
vrs.removeNoGT(f.sample_id)
vrs.removeHomRef(f.sample_id)
vrs.removeFailingFilter()
# create temp vcf file
input_file_bn = os.path.splitext(
    os.path.basename(f.sample_vcf))[0]
sample_vcf = os.path.join(tmp_dir, f.sample_id + '.vcf')
# outp_tsv = os.path.join(tmp_dir, f.sample_id + '.tsv')

# add header lines that are used for de novo analysis
additional_header = vrs.readVcfHeader(
    os.path.join(os.path.dirname(script_name),
                 'header_extra.txt'))
both_head = vrs.vcf_header.split('\n') +\
                           additional_header.split('\n')
both_head = [i for i in both_head if i != '']
vrs.vcf_header = '\n'.join(both_head) + '\n'

vrs.saveAsVcf(sample_vcf)

# annotate with VEP and snpEff
logging.info('preparing to run %s' % script_name)

cmd = ' '.join([script_name,
                sample_vcf,
                os.path.dirname(script_name),
                child_id,
                targ_bed,
                incl_make])
logging.info('Executing \n %s' % cmd)
res = func.runInShell(cmd, True)


