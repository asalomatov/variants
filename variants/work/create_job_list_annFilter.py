import sys
import os

# # batches = 'b1-2 b3 b4 b5 b6 b7 b8 b9'.split()
snp_score = '0.06'
indel_score = '0.0008'
batches = ['b10']
var_types = ['SNP', 'INDEL']
callers = ['hc', 'fb', 'pl']
# callers = ['hc']
work_dir = '/mnt/ceph/users/asalomatov/spark/denovo/def1/'
job_list = []

for b in batches:
    for clr in callers:
        cmd_indel = 'python ~/projects/variants/variants/annotateAndFilter.py %(work_dir)s/indel/%(clr)s/%(b)s/%(clr)s-INDEL-class.csv\
        %(work_dir)s/config/cfg_spark_indel_%(clr)s.yml %(indel_score)s' % locals()
        cmd_snp = 'python ~/projects/variants/variants/annotateAndFilter.py %(work_dir)s/snp/%(clr)s/%(b)s/%(clr)s-SNP-class.csv\
        %(work_dir)s/config/cfg_spark_snp_%(clr)s.yml %(snp_score)s' % locals()
        job_list.append(cmd_indel)
        job_list.append(cmd_snp)
