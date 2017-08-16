#!/mnt/xfs1/home/asalomatov/miniconda2/bin/python
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

arg_parser.add_argument('--remove_tempfiles',
                        nargs='+',
                        type=str,
                        default='Yes',
                        help='Yes/No, the intermidiary files deleted or not')

args = arg_parser.parse_args()
print(args)
input_file = args.input_file[0]
config_file = args.yaml_config_file[0]
prob_cutoff = args.class_probability_threshold[0]
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

dnvo = pandas.read_csv(input_file)
#dnvo.ix[:, 'pred_labels'] = (dnvo['pred_prob'] > prob_cutoff).astype(int)
#dnvo = dnvo[dnvo.pred_labels == 1]
dnvo.reset_index(inplace=True)
if dnvo.empty:
    sys.exit('No mutations at score %s' % prob_cutoff)

tmp_dir = tempfile.mkdtemp()
print(tmp_dir)
input_file_bn = os.path.splitext(os.path.basename(input_file))[0]
outp_tsv = os.path.join(tmp_dir, input_file_bn + '.tsv')
print(outp_tsv)
func.writeTableAsVcf(dnvo, outp_tsv)
# script_name = os.path.basename(os.path.realpath(sys.argv[0]))
script_name = os.path.abspath(pkg_resources.resource_filename(
    'variants',
    'vcf2table.sh'))
#    'vcf2table_notarg.sh'))
cmd = ' '.join([script_name,
                outp_tsv,
                os.path.dirname(script_name),
                input_file_bn,
                targ_bed,
                incl_make,
                vep_refseq])
print(cmd)
func.runInShell(cmd)
vn = summarizeOtherVariants.summarizeMutations(
    os.path.join(tmp_dir,
                 input_file_bn +
                 '-ann.vcf.onePline.tsv'),
    os.path.join(tmp_dir,
                 input_file_bn +
                 '-vep.tsv'),
    input_file_bn,
    output_dir,
    config_file)
if rm_tmp == 'yes':
    cmd = 'rm -rf %s' % tmp_dir
    func.runInShell(cmd)

sys.exit('stop')


def get_spID(x, lab2sp_dict):
    if x[:2] == 'SP':
        return x
    else:
        return lab2sp_dict[x]


labID2spID = func.readYml('/mnt/scratch/asalomatov/data/SPARK/info/baylor_id2sp_descr.yml')
spID2labID = func.readYml('/mnt/scratch/asalomatov/data/SPARK/info/sp_descr2baylor_id.yml')
other_ann = pandas.read_csv('/mnt/xfs1/home/asalomatov/projects/spark/other_calls_b1-4_ALL_ALL_DENOVO.csv')
other_ann['SP_id'] = other_ann.ind_id.apply(lambda z: get_spID(z, labID2spID))
other_ann['lab_id'] = other_ann.SP_id.apply(lambda z: spID2labID[z])
other_ann['v_id'] = other_ann.lab_id.astype(str) + '_' +\
                    other_ann.CHROM.astype(str) + '_' +\
                    other_ann.POS.astype(str)
#other_ann['var_id'] = other_ann.ind_id.astype(str) + '_' +\
#                      other_ann.CHROM.astype(str) + '_' +\
#                      other_ann.POS.astype(str) + '_' +\
#                      other_ann.var_type
other_ann['BCM'] = other_ann.v_id.apply(lambda z: z in bcm_set)
other_ann['Codified'] = other_ann.v_id.apply(lambda z: z in bcm_set)
other_ann['SY'] = other_ann.v_id.apply(lambda z: z in bcm_set)
other_ann['SF'] = False
other_ann['sort_cat'] = 10
other_ann.ix[other_ann.impact_lof == 'True', 'sort_cat'] = 1
other_ann.ix[other_ann.dmg_miss == 'True', 'sort_cat'] = 2
other_ann.ix[(other_ann.coding_var == 'True') & (other_ann.dmg_miss == 'False') &
             (other_ann.missense == 'True'), 'sort_cat'] = 3
other_ann.ix[(other_ann.coding_var == 'True') & (other_ann.sort_cat > 3),
             'sort_cat'] = 4

all_calls_research = pandas.concat([sf_calls_research, other_ann])
all_calls_research = all_calls_research.sort_values(['sort_cat', 'batch', 'lab_id',
                                            'CHROM', 'POS'])
