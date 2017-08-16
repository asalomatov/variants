import sys, os
import pandas
from variants import func
import ped
import argparse
import collections
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants/work/')
from init_baylor_info import spID2labID
from init_baylor_info import labID2spID
from init_baylor_info import lab2batch
from bcm_cdfd_sf_compare import *


other_ann = pandas.read_csv('/mnt/ceph/users/asalomatov/spark/other_calls_b1-2_b3_b4_b5_ALL_ALL_DENOVO.csv')
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
other_ann['Cdfd'] = other_ann.v_id.apply(lambda z: z in codif_set)
other_ann['Clmb'] = other_ann.v_id.apply(lambda z: z in sy_set)
other_ann['SF'] = False
other_ann['batch'] = 'b' + other_ann.lab_id.apply(
    lambda i: str(lab2batch[i]))
other_ann['sort_cat'] = 10
other_ann.ix[other_ann.impact_lof, 'sort_cat'] = 1
other_ann.ix[other_ann.dmg_miss, 'sort_cat'] = 2
other_ann.ix[(other_ann.coding_var) & (other_ann.dmg_miss) &
             (other_ann.missense), 'sort_cat'] = 3
other_ann.ix[(other_ann.coding_var) & (other_ann.sort_cat > 3),
             'sort_cat'] = 4
other_ann['gene'] = other_ann['ANN[*].GENE']
other_ann = other_ann.merge(asd_gene_prob_df, on='gene', how='left')
other_ann = other_ann.merge(ios_anno_df, on='gene', how='left')

all_calls_research = pandas.concat([sf_calls_research, other_ann])
all_calls_research = all_calls_research.sort_values(['sort_cat', 'batch', 'lab_id',
                                        'CHROM', 'POS'])
all_calls_research['N_centers'] = all_calls_research.SF.astype(int) +\
                                  all_calls_research.BCM.astype(int) +\
                                  all_calls_research.Cdfd.astype(int) +\
                                  all_calls_research.Clmb.astype(int)
all_calls_research['var_type'] = all_calls_research.apply(varType, axis=1)


def concatCenters(row):
    res = []
    if row['SF']: res.append('SF')
    if row['BCM']: res.append('BCM')
    if row['Cdfd']: res.append('Cdfd')
    if row['Clmb']: res.append('Clmb')
    return '_'.join(res)

all_calls_research['centers'] = all_calls_research.apply(
    concatCenters, axis=1)
all_calls_research.columns = [i.replace('ANN[*].', '') for i in all_calls_research.columns]
cols_to_output = cols_to_output + ['centers', 'N_centers']
cols_to_output = [i.replace('ANN[*].', '') for i in cols_to_output]
    
