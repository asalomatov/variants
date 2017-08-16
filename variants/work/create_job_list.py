import sys
import os
import glob


denovo_def = sys.argv[1]
kids_files = [i for f in sys.argv[2:] for i in glob.glob(f)]
kids_files.sort()
print 'kids files '
print '\n'.join(kids_files)
#batches = ['b1-2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9']

job_list = []

for kf in kids_files:
    kids = []
    batch_num = os.path.basename(kf).split('_')[0]
    print 'batch %s' % batch_num
    with open(kf, 'r') as f:
        for line in f:
            kids.append(line.strip())
        print 'Number of children %s' % len(kids)
    for k in kids:
        for clr in ['hc', 'fb', 'pl']:
            for var_type in ['snp', 'indel']:
                cmd = 'call_de_novo.py %(k)s /mnt/ceph/users/asalomatov/spark/denovo/def%(denovo_def)s/config/cfg_spark_%(var_type)s_%(clr)s.yml 0.0' % locals()
#                print cmd
                job_list.append(cmd)
print 'Total number of jobs is %s' % len(job_list)
