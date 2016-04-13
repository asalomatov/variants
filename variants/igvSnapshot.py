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
                        nargs='+',
                        help='input file, a tab delimited file containing\
                        ind_id, CHROM, POS columns')
arg_parser.add_argument('ped_file',
                        nargs='+',
                        help='ped file, that can be read by Ped instance')
arg_parser.add_argument('bam_pattern',
                        nargs='+',
                        help='a pattern that describes how to locate bam\
                        file for a ind_id in ped, e.g.,\
                        path/to/%s.bam, where s is ind_id')
arg_parser.add_argument('output_dir',
                        nargs='+',
                        help='igv script and snapshots go to this dir')
arg_parser.add_argument('--igv',
                        type=str,
                        default='/mnt/xfs1/bioinfoCentos7/software/installs/igv/20151202/igv.jar',
                        help='path to igv.jar to overwrite\
                        the hardcoded default')

args = arg_parser.parse_args()
print args
output_dir = args.output_dir[0]

myped = ped.Ped(args.ped_file[0])
myped.addBam(field='ind_id', file_pat=args.bam_pattern[0])
myped.ped.dropna(subset=['bam'], inplace=True)
myped.ped.reset_index(inplace=True)
myped.addBai(field='ind_id', file_pat=args.bam_pattern[0] + '.bai')
myped.ped.dropna(subset=['bai'], inplace=True)
myped.ped.reset_index(inplace=True)

sample_bam = 'bambam'
chr_pos = '1:100'
sample_snapshot_name = 'chacha.png'

tmpl1 = """

igv
new
genome hg19
snapshotDirectory %(output_dir)s
load %(sample_bam)s
goto %(chr_pos)s
snapshot %(sample_snapshot_name)s

"""

tmpl3 = """

igv
new
genome hg19
snapshotDirectory %(output_dir)s
load %(sample_bam)s
load %(father_bam)s
load %(mother_bam)s
goto %(chr_pos)s
snapshot %(sample_snapshot_name)s

"""



