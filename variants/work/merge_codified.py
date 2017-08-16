import os
import sys
import pandas
import glob

inp_dir = sys.argv[1:]
d_l = []
for i_d in inp_dir:
    print 'dir is ', i_d
    f_l = glob.glob(os.path.join(i_d, '*.tsv'))
    print 'file list is ', f_l
    for f in f_l:
        try:
            d = pandas.read_table(f, dtype=str)
        except:
            continue
        d['ind_id'] = os.path.basename(f).split('.')[0]
        d_l.append(d)


df = pandas.concat(d_l,)
df.reset_index(inplace=True, drop=True)


def extractAlleles(x):
    try:
        y = x.split('>')
    except:
        return pandas.Series([None, None], ['REF', 'ALT'])
    return pandas.Series([y[0][-1], y[1][0]], ['REF', 'ALT'])

def extractChrPos(x):
    y = x.split(':')
    return pandas.Series([y[0], y[1]], ['CHROM', 'POS'])


df = df.merge(df.Position.apply(extractChrPos), left_index=True, right_index=True)
#df = df.merge(df['Transcript data'].apply(extractAlleles), left_index=True, right_index=True)
df['REF'] = None
df['ALT'] = None
df = df[['CHROM', 'POS', 'REF', 'ALT', 'Mutation Type', 'Gene', 'ind_id']]
df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'effect', 'gene', 'ind_id']
df['var_id'] = df.ind_id.astype(str) + '_' + df.CHROM.astype(str) + '_' + df.POS.astype(str)
#d.columns = ['gene', 'pos', 'REF', 'ALT', 'effect', 'gene']
