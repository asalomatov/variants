
import os
import variants
import pysam
import subprocess
import sys
import tempfile
import errno
import glob

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

def listFiles(pattern):
    l = glob.glob(pattern)
    if len(l) == 0:
        return None
    elif len(l) == 1:
        return l[0]
    else:
        return l
    
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

def runInShell(cmd, return_output=False):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode == 0:
        if return_output:
            return stdout
        else:
            return 0
    else:
        sys.exit(stderr)
       
  
def runBamReadcounts(vcffile, bamfile, output_dir, genome_ref = \
    '/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'):
    """exctract regions from vcf and supply them to bam-readcount"""
    func.makeDir(output_dir)
    if func.checkFile(vcffile) is None:
        sys.exit('No vcf file!')
    if func.checkFile(bamfile) is None:
        sys.exit('No vcf file!')
    cat = 'cat '
    if os.path.splitext(vcffile)[1] == '.gz':
        cat = 'zcat '
    bam_name = os.path.splitext(os.path.basename(bamfile))[0]
    tmp_bed = tempfile.mktemp(suffix='.bed', prefix=bam_name+'_')    
    outp_fn = os.path.join(output_dir, bam_name+'.txt')    
    cmd1 = cat + vcffile + \
        " | grep -v ^# | awk \'{print $1\"\t\"$2\"\t\"$2}' > " + tmp_bed
    cmd2 = ' '.join(['bam-readcount -f', genome_ref, bamfile, '-l',tmp_bed, \
            '2>/dev/null 1>'+outp_fn])
    cmd3 = 'rm '+tmp_bed
    cmd = ';'.join([cmd1, cmd2, cmd3])
    return func.runInShell(cmd)
    
def runBamReadcountsRegions(regionsfile, bamfile, output_dir, genome_ref = \
    '/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'):
    """exctract regions from vcf and supply them to bam-readcount"""
    func.makeDir(output_dir)
    if func.checkFile(regionsfile) is None:
        sys.exit('No vcf file!')
    if func.checkFile(bamfile) is None:
        sys.exit('No vcf file!')
    cat = 'cat '
    bam_name = os.path.splitext(os.path.basename(bamfile))[0]
    outp_fn = os.path.join(output_dir, bam_name+'.txt')    
    cmd = ' '.join(['bam-readcount -f', genome_ref, bamfile, '-l',regionsfile, \
            '2>/dev/null 1>'+outp_fn])
    return func.runInShell(cmd)
    
 
    
    
    
    
 