#!/bin/env python

import os
import sys
import pandas
import glob

inp_dir = sys.argv[1:]
d_l = []
for i_d in inp_dir:
    print 'dir is ', i_d
    f_l = glob.glob(os.path.join(i_d, '*.CNV.report'))
    print 'file list is ', f_l
    for f in f_l:
        try:
            d = pandas.read_table(f, dtype=str)
        except:
            continue
#        d.columns = ['CHROM', 'POS', 'REF', 'ALT', 'effect', 'gene']
        d['ind_id'] = os.path.basename(f).split('.')[0]
        d_l.append(d)

df = pandas.concat(d_l,)
df.reset_index(inplace=True, drop=True)
