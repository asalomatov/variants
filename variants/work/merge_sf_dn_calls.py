#!/usr/bin/env python
import sys
import pandas

dn_calls = sys.argv[1:]
dn_list = []
for i in dn_calls:
    x = pandas.read_csv(i)
    x['CHROM'] = x.CHROM.astype(str)
    x['POS'] = x.POS.astype(int)
    if 'hc-' in i and 'jhc-' not in i:
        x['clr'] = 'hc'
    elif 'jhc-' in i:
        x['clr'] = 'jhc'
    elif 'fb-' in i:
        x['clr'] = 'fb'
    elif 'pl-' in i:
        x['clr'] = 'pl'
    else:
        x['clr'] = 'ukn'
    x['var_id'] = x.ind_id + '__' +\
                  x.CHROM + '__' +\
                  x.POS.astype(str)
    dn_list.append(x)
dn = dn_list[0]
#cols_to_merge_on = list(dn.columns[:-1])
cols_to_merge_on = ['var_id']
orig_dtypes = dn[cols_to_merge_on].dtypes.to_dict()
print orig_dtypes
for i in dn_list[1:]:
    print 'merging'
    dn = dn[cols_to_merge_on + ['clr']].merge(i[cols_to_merge_on + ['clr']],
                                              on=cols_to_merge_on,
                                              how='outer')
    for i in orig_dtypes:
        dn[i] = dn[i].astype(orig_dtypes[i])
