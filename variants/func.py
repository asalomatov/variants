
import os
import variants
import pysam

def makeDir(path):
    try:
        os.mkdir(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def checkFile(fname):
    """file name if file exists, None otherwise"""
    if os.path.isfile(fname):
        return fname
    else:
        return None

def run_once(f):
    def wrapper(*args, **kwargs):
        if not wrapper.has_run:
            wrapper.has_run = True
            return f(*args, **kwargs)
    wrapper.has_run = False
    return wrapper

def fileHasVariant(fname, fam_id, chrom, pos_start, pos_end):
    myvars = variants.Variants(fname, fam_id, chrom, pos_start, pos_end)
    next_var = None
    fam_var = []
    for next_var in myvars.vcf_reader:
        alt_allel = []
        for nucl_alt in next_var.ALT: 
            alt_allel.append(nucl_alt.sequence)
        for v in alt_allel:
            fam_var.append('_'.join([fam_id, next_var.CHROM, str(next_var.POS), next_var.REF, v]))
    return fam_var

def refAtPos(chrom, pos, genref='/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'):
    ref_allel = pysam.faidx(genref, str(chrom)+':'+str(pos)+'-'+str(pos))[1].strip()
    return ref_allel

def df2sklearn(mydf, col_to_keep):
    if 'status' in mydf.columns:
        mydf['status01'] = 1
        mydf['status01'][mydf['status'] == 'N'] = 0
        col_to_keep += ['status01']
    col_to_keep = list(set(col_to_keep).intersection(set(mydf.columns)))
    print col_to_keep
    #res = mydf[col_to_keep]
    mydf[col_to_keep] = mydf[col_to_keep].astype(float)
    mydf = mydf.dropna(subset = col_to_keep)
    return mydf[col_to_keep] 

def varType(row):
    if len(row['Ref']) == len(row['Alt']):
        return 'snp'
    else:
        return 'indel'

