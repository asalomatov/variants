import sys
#sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import ped
import func
import pandas
import numpy
import os
import argparse

arg_parser = argparse.ArgumentParser(
    description='Create a batch script for igv,\
                run igv in batch mode to generate\
                snapshots')
arg_parser.add_argument('input_file',
                        help='input file, a tab delimited file containing\
                        ind_id, CHROM, POS , descr columns')
arg_parser.add_argument('ped_file',
                        help='ped file, that can be read by Ped instance,\
                        should have bam column following the mandotary cols')
arg_parser.add_argument('output_dir',
                        help='igv script and snapshots go to this dir')
arg_parser.add_argument('--igv',
                        type=str,
                        default='java -jar /mnt/xfs1/bioinfoCentos7/software/installs/igv/20151202/igv.jar',
                        help='path to igv.jar to overwrite\
                        the hardcoded default')
arg_parser.add_argument('--genome_build',
                        type=int,
                        default=19,
                        help='19, 37, or 38')

args = arg_parser.parse_args()
print args
input_file = args.input_file
output_dir = args.output_dir
if args.genome_build in [19, 37]:
    genome_b = 'hg19'
elif args.genome_build in [38]:
    genome_b = 'genome.fa'
else:
    sys.exit('Unknown genome build, use 19, 37, or 38')
func.runInShell('mkdir -p ' + output_dir)

igv_inp = pandas.read_csv(input_file, dtype=str)
#igv_inp = pandas.read_table(input_file, dtype=str)
#for reg
if args.genome_build in [38]:
    igv_inp.SP_id = igv_inp.ind_id
    igv_inp.lab_id = igv_inp.ind_id

myped = ped.Ped(args.ped_file, ['bam', 'vcf'])
#myped = ped.Ped(args.ped_file, ['batch', 'bam'])


tmpl1 = """

new
genome %(genome_b)s
snapshotDirectory %(output_dir)s
load %(sample_bam)s
goto %(chr_pos)s
snapshot %(sample_snapshot_name)s

"""

tmpl3 = """

new
genome %(genome_b)s
snapshotDirectory %(output_dir)s
load %(sample_bam)s
load %(father_bam)s
load %(mother_bam)s
goto %(chr_pos)s
snapshot %(sample_snapshot_name)s

"""

snapshot_list = []

for i, row in igv_inp.iterrows():
    smpl_id = row['SP_id']
    trio_id = myped.getFamily(smpl_id)
    print trio_id
    sample_bam = myped.getIndivBAM(smpl_id)
    father_bam = myped.getFaBam(trio_id)
    mother_bam = myped.getMoBam(trio_id)
    chr_pos = ':'.join([row['CHROM'],
                        '-'.join([str(int(row['POS']) - 50),
                                  str(int(row['POS']) + 50)])])
    sample_snapshot_name = '_'.join([row['SP_id'],
                                     row['lab_id'],
                                     row['CHROM'],
                                     row['POS'],
                                     row['GENE'],
                                     row['centers']]) + '.png'
    print sample_bam, father_bam, mother_bam, chr_pos, sample_snapshot_name
    snapshot_list.append(tmpl3 % locals())

igv_scr = '\n'.join(snapshot_list)

func.runInShell('mkdir -p ' + output_dir)
input_file_bn = os.path.basename(input_file)
script_out = os.path.join(output_dir, input_file_bn + '.igv')
with open(script_out, 'w') as f:
    f.write(igv_scr)

cmd = ' '.join([args.igv, '-b', script_out])
print cmd
#func.runInShell(cmd)

sys.exit(1)


