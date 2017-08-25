from __future__ import print_function
import sys
from variants import func
import pandas
import os
import yaml
import summarizeOtherVariants
import argparse
import pkg_resources
import tempfile

print('command line:')
print(sys.argv)

arg_parser = argparse.ArgumentParser(
    description='Annotate and filter de novo mutations')
arg_parser.add_argument('input_vcf',
                        nargs='+',
                        type=str,
                        help='Annotate a vcf')
arg_parser.add_argument('yaml_config_file',
                        nargs='+',
                        type=str,
                        help='A config file, defining necessary variables.')
arg_parser.add_argument('--remove_tempfiles',
                        nargs='+',
                        type=str,
                        default='Yes',
                        help='Yes/No, the intermidiary files deleted or not')

args = arg_parser.parse_args()
print(args)
input_file = args.input_vcf[0]
config_file = args.yaml_config_file[0]
rm_tmp = args.remove_tempfiles[0].lower()


# get parameters from yaml config file
with open(config_file, 'r') as f:
    cfg = yaml.safe_load(f)
min_DP = cfg['min_DP']
var_type = cfg['variant_type']
vcf_pat = cfg['vcf_pattern']
ped_file = cfg['ped_file']
ped_file_extended = cfg['ped_file_extended']
known_vars = cfg['known_variants']
output_dir = cfg['output_directory']
targ_bed = cfg['target_bed']
genome_build = int(cfg['genome_build'])
vep_refseq =  cfg['vep_refseq']
if genome_build == 19 or genome_build == 37:
    incl_make = '/mnt/xfs1/home/asalomatov/projects/pipeline/ppln/include.mk'
elif int(genome_build) == 38:
    incl_make = '/mnt/xfs1/home/asalomatov/projects/pipeline/ppln/include_hg38.mk'
else:
    sys.exit('Only builds 19, 37, 38 are supported')

tmp_dir = tempfile.mkdtemp()
print(tmp_dir)
input_file_bn = os.path.splitext(os.path.basename(input_file))[0]
input_lile_dir = os.path.dirname(input_file)
script_name = os.path.abspath(pkg_resources.resource_filename(
    'variants',
    'vcf2tablee.sh'))
script_dir = os.path.dirname(script_name)

cmd = """
vcfintersect -b %(targ_bed)s %(input_file)s > %(tmp_dir)s/%(input_file_bn)s.vcf
echo 'running VEP'
make -f %(script_dir)s/annVEP.mk INCLMK=%(incl_make)s VEPREFSEQ=%(vep_refseq)s PREFIX=%(input_file_bn)s SUFFIX=.vcf INDIR=%(tmp_dir)s OUTDIR=%(tmp_dir)s
"""

print(cmd % locals())

func.runInShell(cmd % locals())

if rm_tmp == 'yes':
    cmd = 'rm -rf %s' % tmp_dir
    func.runInShell(cmd)

