import os
#import pysam
import subprocess
import sys
import tempfile
import errno
import glob
import pandas
import numpy
import collections
#import train
import pysam
import yaml


vcf_required_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
                                'QUAL', 'FILTER']

def mergeFieldsForVariant(x):
    y = (x.HGVSc + ';' + x.EXON + ';' + x.INTRON + ';' + x.HGVSp).tolist()
    y = [i for i in y if i[0] != '-']
    return '|'.join(y)



def parentForOffspring(smpl_id, ped_file, sex=1):
    ped = pandas.read_table(ped_file, usecols=range(6), header=None, dtype=str)
    ped.columns = ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                   'sex', 'pheno']
    if sex == 1:
        fa_id = ped.fa_id[ped.ind_id == smpl_id].iloc[0]
        return fa_id
    else:
        mo_id = ped.mo_id[ped.ind_id == smpl_id].iloc[0]
        return mo_id


def siblingsForOffspring(smpl_id, ped_file):
    ped = pandas.read_table(ped_file, usecols=range(6), header=None, dtype=str)
    ped.columns = ['fam_id', 'ind_id', 'fa_id', 'mo_id',
                   'sex', 'pheno']
    ped['ff'] = ped.fa_id + '_' + ped.mo_id
    ff = ped.ff[ped.ind_id == smpl_id].iloc[0]
    # print(ped.ff.value_counts().head())
    prnts_prob = ff.split('_') + [smpl_id]
    children = ped.ind_id[ped.ff == ff].tolist()
    sibs = [i for i in children if i not in prnts_prob]
    return '_'.join(sibs)


def sibFromDataMart(smpl_id, dmrt_df):
	sibs = dmrt_df.ix[dmrt_df.user_sp_id == smpl_id, 'sibling_sp_id'].tolist()
	sibs_cons = dmrt_df.ix[(dmrt_df.user_sp_id == smpl_id) &
		(dmrt_df.sibling_genetic_consented == 1), 'sibling_sp_id'].tolist()
	sibs_dna = dmrt_df.ix[(dmrt_df.user_sp_id == smpl_id) &
		(dmrt_df.sibling_with_dna == 1), 'sibling_sp_id'].tolist()
	sibs_seq = dmrt_df.ix[(dmrt_df.user_sp_id == smpl_id) &
		(dmrt_df.sibling_sequenced == 1), 'sibling_sp_id'].tolist()
	sibs_aff = dmrt_df.ix[(dmrt_df.user_sp_id == smpl_id) &
		(dmrt_df.sibling_affected == 1), 'sibling_sp_id'].tolist()
	sibs_unaff = [i for i in sibs if i not in sibs_aff]
	twins = dmrt_df.ix[(dmrt_df.user_sp_id == smpl_id) &
		(dmrt_df.twins == 1), 'sibling_sp_id'].tolist()
	twins_unaff = [i for i in twins if i in sibs_unaff]
	twins_aff = [i for i in twins if i in sibs_aff]
	res = pandas.Series(
					['_'.join(sibs),
                    '_'.join(sibs_cons),
                    '_'.join(sibs_dna),
                    '_'.join(sibs_seq),
                    '_'.join(sibs_aff),
					'_'.join(sibs_unaff),
					'_'.join(twins_aff),
					'_'.join(twins_unaff)],
					['Siblings',
                    'Consented siblings',
                    'Siblings with dna',
                    'Sequenced siblings',
                    'Affected siblings',
					'Unaffected siblings',
					'Affected twins',
					'Unaffected twins'])
	return res


def sumGene(x):
    Count = len(x)
    Count_LOF = sum(x.lof.astype(str) == 'True')
    Count_MIS = sum(x.missense.astype(str) == 'True')
    N_indiv = round(len(x.SP_id.unique()), 0)
    indiv = '_'.join(x.SP_id.unique().tolist())
    SFARI_score = x.SFARIscore.iloc[0]
    lof_z_perc_rank = round(x.lof_z_perc_rank.iloc[0], 3)
    mis_z_perc_rank = round(x.mis_z_perc_rank.iloc[0], 3)
    pLI_perc_rank = round(x.pLI_perc_rank.iloc[0], 3)
    asd_score_perc_rank = round(x.asd_score_perc_rank.iloc[0], 3)
    LGDscore_perc_rank = round(x.LGDscore_perc_rank.iloc[0], 3)
    RVIS_perc_rank = round(x.RVIS_perc_rank.iloc[0], 3)
    return pandas.Series(
        [Count, Count_LOF, Count_MIS, N_indiv,
         SFARI_score,
         lof_z_perc_rank, mis_z_perc_rank,
         pLI_perc_rank,
         asd_score_perc_rank,
         LGDscore_perc_rank,
         RVIS_perc_rank,
         indiv],
        index=['N', 'N_LOF', 'N_DMIS',
               'N_indiv', 'SFARI',
               'lof_z', 'mis_z',
               'pLI',
               'asd(Olga)',
               'LGD',
               'RVIS',
               'indiv']
    )


def sumGeneGPF(x):
    Count = len(x)
    Count_LGD = sum(x['worst requested effect'] != 'missense')
    Count_MIS = sum(x['worst requested effect'] == 'missense')
    return pandas.Series(
        [Count, Count_LGD, Count_MIS],
        index=['Number of mutations', 'Number of LGD', 'Number of MIS']
    )


# def sumGene(x):
#     Count = len(x)
#     Count_LOF = sum(x.missense.astype(str) == 'True')
#     Count_MIS = sum(x.lof.astype(str) == 'True')
#     N_indiv = round(len(x.SP_id.unique()), 0)
#     indiv = '_'.join(x.SP_id.unique().tolist())
#     SFARI_score = x.SFARIscore.iloc[0]
#     lof_z_perc_rank = round(x.lof_z_perc_rank.iloc[0], 3)
#     mis_z_perc_rank = round(x.mis_z_perc_rank.iloc[0], 3)
#     pLI_perc_rank = round(x.pLI_perc_rank.iloc[0], 3)
#     asd_score_perc_rank = round(x.asd_score_perc_rank.iloc[0], 3)
#     LGDscore_perc_rank = round(x.LGDscore_perc_rank.iloc[0], 3)
#     RVIS_perc_rank = round(x.RVIS_perc_rank.iloc[0], 3)
#     return pandas.Series(
#         [Count, N_indiv,
#          SFARI_score,
#          lof_z_perc_rank, mis_z_perc_rank,
#          pLI_perc_rank,
#          asd_score_perc_rank,
#          LGDscore_perc_rank,
#          RVIS_perc_rank,
#          indiv],
#         index=['Count', 'N_indiv',
#                'SFARI_score',
#                'lof_z_perc_rank', 'mis_z_perc_rank',
#                'pLI_perc_rank',
#                'asd_score_perc_rank',
#                'LGDscore_perc_rank',
#                'RVIS_perc_rank',
#                'indiv']
#     )


def numCoding(df, p_snp, p_indel):
    snp = df[(df.pred_prob.astype(float) > p_snp) &
               (df.VARTYPE == 'SNP')]
    indel = df[(df.pred_prob.astype(float) > p_indel) &
               (df.VARTYPE != 'SNP') &
               ((df.inherit_prnts.astype(int) == 0) |
                (df.alt_DP_fa.astype(int) + df.alt_DP_mo.astype(int) < 3))]
    print('num of coding snp: %s' % snp[snp.effect_cat != 'other'].shape[0])
    print(snp.inherit_prnts.astype(int).value_counts())
    print('num of coding indel: %s' %indel[indel.effect_cat != 'other'].shape[0])
    print(indel.inherit_prnts.astype(int).value_counts())
    return pandas.concat([snp, indel])


def varId(df, idfield='SP_id'):
    varid = df[idfield].astype(str) + '_' +\
            df.CHROM.astype(str) + '_' +\
            df.POS.astype(str)
    return varid


def checkConcordance(call_set1, call_set2, idfield='SP_id'):
    df1 = pandas.read_csv(call_set1)
    df2 = pandas.read_csv(call_set2)
    set1 = set(varId(df1[df1.effect_cat != 'other'], idfield))
    set2 = set(varId(df2[df2.effect_cat != 'other'], idfield))
    both = len(set1.intersection(set2))
    set1not2 = len(set1.difference(set2))
    set2not1 = len(set2.difference(set1))
    print(call_set1 + ' only: %s' % set1not2)
    if len(set1.difference(set2)) < 6:
        print(set1.difference(set2))
    print('both sets: %s' % both)
    print(call_set2 + ' only: %s' % set2not1)
    if len(set2.difference(set1)) < 6:
        print(set2.difference(set1))


def addRank(df, clm_name):
    df[clm_name + '_rank'] = df[clm_name].rank()


def addPercRank(df, clm_name):
    df[clm_name + '_perc_rank'] = df[clm_name].rank(pct=True)


def readVcfToDF1(fname, sample_list=None, chunk_size=None):
    """read vcf file into pandas DF without parsing.
    If sample_list is None, it'll read all of the samples"""
    # after github/hammerlab/varcode/vcf but keeping sample information
    path = fname
    compression = None
    if path.endswith(".gz"):
        compression = "gzip"
    elif path.endswith(".bz2"):
        compression = "bz2"
    cat = 'cat'
    if compression is not None:
        cat = 'zcat'
    cmd = ' '.join([cat, path, '| head -10000 | grep ^# | grep -v ^##'])
    vcf_clmns = runInShell(cmd, True).split('\t')
    vcf_clmns = [x.strip() for x in vcf_clmns]
    vcf_clmns = [x.strip('#') for x in vcf_clmns]
    print(vcf_clmns)
    df_cols = []
    df_cols_ind = []
    if sample_list is None:
        df_cols = vcf_clmns
        df_cols_ind = range(len(df_cols))
    else:
        df_cols = vcf_required_fields[:]
        df_cols_ind = range(len(df_cols))
        smp_indexes = []
        for smp in sample_list:
            smp_ind = vcf_clmns.index(smp)
            smp_indexes.append(smp_ind)
        smp_indexes.sort()
        print('smp indexes:')
        print smp_indexes
        for i in smp_indexes:
            print('appending')
            print(vcf_clmns[i])
            df_cols.append(vcf_clmns[i])
            df_cols_ind.append(i)
    df_field_types = collections.OrderedDict()
    for i in df_cols:
        df_field_types[i] = str
    df_field_types['POS'] = int
    print('reading columns')
    print(df_cols)
    reader = pandas.read_table(
        path,
        compression=compression,
        comment="#",
        chunksize=chunk_size,
        dtype=df_field_types,
        names=df_cols,
        usecols=df_cols_ind)
    return reader


def readVcfToDF(fname, chunk_size=None):
    """read vcf or file into pandas DF without parsing,
    also works for reading VEP output files"""
    # after github/hammerlab/varcode/vcf but keeping sample information
    compression = None
    if fname.endswith(".gz"):
        compression = "gzip"
    elif fname.endswith(".bz2"):
        compression = "bz2"
    cat = 'cat'
    if compression is not None:
        cat = 'zcat'
    cmd = ' '.join([cat, fname, '| grep ^# | grep -v ^##'])
    vcf_clmns = runInShell(cmd, True).split('\t')
    vcf_clmns = [x.strip() for x in vcf_clmns]
    vcf_clmns = [x.strip('#') for x in vcf_clmns]
    vcf_field_types = collections.OrderedDict()
    for i in vcf_clmns:
        vcf_field_types[i] = str
    reader = pandas.read_table(
        fname,
        compression=compression,
        comment="#",
        chunksize=chunk_size,
        dtype=vcf_field_types,
        names=list(vcf_field_types),
        usecols=range(len(vcf_field_types)))
    return reader


def readYml(path):
    with open(path, 'r') as f:
        res = yaml.safe_load(f)
    return res


def dumpYml(path, x):
    with open(path, 'w') as f:
        yaml.dump(x, f, default_flow_style=False)


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


def refAtPos(chrom, pos, genref):
    ref_allel = pysam.faidx(genref,
                            str(chrom)+':'+str(pos)+'-'+str(pos))[1].strip()
    return ref_allel

# /mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa

def bamrcIndel2vcfIndel(ind_id, chrom, position, ref, smpl_alt, genref):
    REF = ref
    ALT = smpl_alt.split('_')[0]
    POS = position
    c_ins = ALT[0] == '+'
    c_del = ALT[0] == '-'
    ALT = ALT.strip('+')
    ALT = ALT.strip('-')
    if c_ins:
        ALT = REF + ALT
    elif c_del:
        REF = ALT
        # POS -= 1  commented out for train set use only
        ref_at_pos = refAtPos(chrom,
                              POS,
                              genref)
        REF = ref_at_pos + REF
        ALT = ref_at_pos
    else:
        sys.exit('Unknown type')
    res = pandas.Series([ind_id, chrom, POS, REF, ALT],
                        ['ind_id', 'CHROM', 'POS', 'REF', 'ALT'])
    return res

def vepVar2vcfVar(vep_tsv_row, genref):
    """
    Convert VEP variant description (Uploaded_variation),
    e.g., 3_148802449_-/TTTAG,
    to VCF fields, CHROM, POS, REF, ALT.
    """
    # print(vep_tsv_row.Uploaded_variation)
    if vep_tsv_row.Uploaded_variation[:2] != 'rs':
        # print('no rs id')
        chrom, pos, refalt = vep_tsv_row.Uploaded_variation.split('_')
        ref, alt = refalt.split('/')
        if ref != '-' and alt != '-':
            #  snp
            pass
        elif ref == '-' and alt != '-':
            #  insertion
            pos = int(pos) - 1
            ref = refAtPos(chrom=chrom, pos=pos, genref=genref)
            alt = ref + alt
        elif ref != '-' and alt == '-':
            #  deletion
            pos = int(pos) - 1
            alt = refAtPos(chrom=chrom, pos=pos, genref=genref)
            ref = alt + ref
        else:
            sys.exit('Unknown mutation type: %s' % vep_var)
    else:
        chrom, pos = vep_tsv_row.Location.split(':')
        if len(pos.split('-')) > 1:
            sys.exit('this is not SNP!')
        alt = vep_tsv_row.Allele
        ref = refAtPos(chrom=chrom, pos=pos, genref=genref)
    res = pandas.Series([chrom, pos, ref, alt],
                        ['CHROM', 'POS', 'REF', 'ALT'])
    return res


def df2sklearn(mydf, col_to_keep):
    if 'status' in mydf.columns:
        mydf['status01'] = 1
        mydf['status01'][mydf['status'] == 'N'] = 0
        col_to_keep += ['status01']
    col_to_keep = list(set(col_to_keep).intersection(set(mydf.columns)))
    print col_to_keep
    #res = mydf[col_to_keep]
    mydf[col_to_keep] = mydf[col_to_keep].astype(float)
    mydf = mydf.dropna(subset=col_to_keep)
    return mydf[col_to_keep]


def varType(row):
    c1 = len(row['Ref']) == len(row['Alt'])
#    c2 = ',' in row['Ref']
#    c3 = ',' in row['Alt']
    if c1:
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
        sys.stderr.write(stderr)
        return 1


def createTempFile(suffix, prefix, tmp_dir=None):
    if not tmp_dir:
        tmp_dir = tempfile.mkdtemp()
    makeDir(tmp_dir)
    tmp_file = tempfile.mktemp(dir=tmp_dir,
                               suffix=suffix,
                               prefix=prefix)
    return (tmp_dir, tmp_file)


def runBamReadcounts(vcffile, bamfile, output_dir,
                     genome_ref='/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/\
                     genomes/Hsapiens/GRCh37/seq/GRCh37.fa',
                     brc='/mnt/scratch/asalomatov/software/installs/\
                     bin/bam-readcount'):
    """exctract regions from vcf and supply them to bam-readcount"""
    makeDir(output_dir)
    if checkFile(vcffile) is None:
        sys.exit('No vcf file!')
    if checkFile(bamfile) is None:
        sys.exit('No bam file!')
    cat = 'cat '
    if os.path.splitext(vcffile)[1] == '.gz':
        cat = 'zcat '
    bam_name = os.path.splitext(os.path.basename(bamfile))[0]
    tmp_dir = tempfile.mkdtemp()
    tmp_bed = tempfile.mktemp(dir=tmp_dir, suffix='.bed', prefix=bam_name+'_')
    outp_fn = os.path.join(output_dir, bam_name+'.txt')
    cmd1 = cat + vcffile + \
        " | grep -v ^# | awk \'{print $1\"\t\"$2\"\t\"$2}' > " + tmp_bed
    cmd2 = ' '.join([brc, '-f', genome_ref, bamfile, '-l',tmp_bed, \
            ' > '+outp_fn])
    cmd3 = ' '  # 'rm -rf '+tmp_dir
    cmd = ';'.join([cmd1, cmd2, cmd3])
    return runInShell(cmd)


def runBamReadcountsRegions(regionsfile,
                            bamfile,
                            output_file,
                            genome_ref='/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa',
                            brc='/mnt/scratch/asalomatov/software/installs/bin/bam-readcount'):
    """exctract features using bam-readcount"""
    if checkFile(regionsfile) is None:
        sys.exit('No regions file!')
    if checkFile(bamfile) is None:
        sys.exit('No bam file!')
    cmd = ' '.join([brc, '-f', genome_ref, bamfile, '-l', regionsfile,
                    ' > ' + output_file])
    return runInShell(cmd)


def readBamReadcount1(file_name, vartype):
    """Read bam-readcount output into pandas DF,
    vartype = ['SNP', 'INS', 'DEL', 'INDEL']
    """
    dir_name = os.path.dirname(file_name)
    tmp_file = file_name + '.' + vartype
    if vartype == 'SNP':
        cmd1 = 'cat ' + file_name + \
               " | grep -v + | grep -v - >" + tmp_file
    elif vartype == 'INS':
        cmd1 = 'cat ' + file_name + \
               " | grep + >" + tmp_file
    elif vartype == 'DEL':
        cmd1 = 'cat ' + file_name + \
               " | grep - >" + tmp_file
    elif vartype == 'INDEL':
        cmd1 = 'cat ' + file_name + \
               " | grep -e '-\|+' >" + tmp_file
    else:
        sys.exit('unknown vartype...')
    runInShell(cmd1)
    if vartype == 'SNP':
        clm_names = ['CHROM', 'POS',
                     'REF', 'DP', 'A', 'C', 'G', 'T']
        clm_dtypes = collections.OrderedDict()
        for i in clm_names:
            clm_dtypes[i] = str
        reader = pandas.read_csv(
            tmp_file,
            header=None,
            dtype=clm_dtypes,
            names=list(clm_names),
            usecols=[0,1,2,3,5,6,7,8],
            sep='\t')
        res = reader.apply(parseBamReadcount, axis=1)
        return reader[['CHROM', 'POS', 'REF', 'DP']].merge(res,
                                                           left_index=True,
                                                           right_index=True)


def readBamReadcount(file_name, vartype='snp', n_clmns_per_allele=14):
    """Read bam-readcount output into pandas DF.
    """
    clm_names = ['CHROM', 'POS', 'REF', 'DP', 'ZERO', 'A', 'C', 'G', 'T',
                 'N', 'INDEL', 'INDEL1', 'INDEL2', 'INDEL3'] + list('12345')
    clm_dtypes = collections.OrderedDict()
    for i in clm_names:
        clm_dtypes[i] = str
    reader = pandas.read_csv(
        file_name,
        header=None,
        dtype=clm_dtypes,
        names=list(clm_names),
#        usecols=range(len(clm_names)),
        sep='\t')
    reader['INDEL'][reader.INDEL.isnull()] = ':'.join(['0'] *
                                                      n_clmns_per_allele)
    reader['INDEL1'][reader.INDEL1.isnull()] = ':'.join(['0'] *
                                                        n_clmns_per_allele)
    reader['INDEL2'][reader.INDEL2.isnull()] = ':'.join(['0'] *
                                                        n_clmns_per_allele)
    reader['INDEL3'][reader.INDEL3.isnull()] = ':'.join(['0'] *
                                                        n_clmns_per_allele)
    reader.drop('ZERO', axis=1, inplace=True)
    reader.drop(list('12345'), axis=1, inplace=True)
    # return reader
    if vartype.lower() == 'snp':
        print('parsing %s' % vartype)
        res = reader.apply(parseBamReadcountSNP, axis=1)
    elif vartype.lower() == 'indel':
        print('parsing %s' % vartype)
        res = reader.apply(parseBamReadcountIndel, axis=1)
    else:
        sys.exit('vartype unrecognized, must be SNP or INDEL')
    return reader[['CHROM', 'POS', 'REF', 'DP']].merge(res, left_index=True,
                                                       right_index=True)


def parseBamReadcountIndel(row):
    """Parsing indels only"""
    # first split ref data
    ref_split = row[row['REF']].split(':')
    ref_split = ref_split[:1] + [float(x) for x in ref_split[1:]]
    # now arrgegate information for alt allels
    possib_alt = ['INDEL', 'INDEL1', 'INDEL2', 'INDEL3']
    num_allels = 0
    alt_read_count = 0
    cum_array = numpy.zeros(12)
    alts = []
    for i in possib_alt:
        row_list = row[i].split(':')
        if int(row_list[1]) > 0:
            num_allels += 1
            allel_count = float(row_list[1])
            alt_read_count += allel_count
            alts += row_list[:2]
            cum_array += allel_count * numpy.array(row_list[2:], dtype=float)
    if alt_read_count > 0:
        cum_array = cum_array / alt_read_count
    clmns_detail = ['base',
                    'count',
                    'avg_mapping_quality',
                    'avg_base_quality',
                    'avg_se_mapping_quality',
                    'num_plus_strand',
                    'num_minus_strand',
                    'avg_pos_as_fraction',
                    'avg_num_mismatches_as_fraction',
                    'avg_sum_mismatch_qualities',
                    'num_q2_containing_reads',
                    'avg_dist_to_q2_start_in_q2_reads',
                    'avg_clipped_length',
                    'avg_dist_to_effective_3p_end']
    res = pandas.Series(ref_split + [num_allels, '_'.join(alts), alt_read_count] +
                        cum_array.tolist(),
                        ['REF_'+ x for x in clmns_detail] + ['num_allels'] +
                        ['ALT_'+ x for x in clmns_detail])
    return res


def parseBamReadcountSNP(row):
    """Parsing SNPs only"""
    # first split ref data
    ref_split = row[row['REF']].split(':')
    ref_split = ref_split[:1] + [float(x) for x in ref_split[1:]]
    #  split indel data
    indel_split = row['INDEL'].split(':')
    indel_split = indel_split[:1] + [float(x) for x in indel_split[1:]]
    # now arrgegate information for alt allels
    possib_alt = []
    if row['REF'] == 'A': possib_alt = ['C', 'G', 'T']
    if row['REF'] == 'C': possib_alt = ['A', 'G', 'T']
    if row['REF'] == 'G': possib_alt = ['A', 'C', 'T']
    if row['REF'] == 'T': possib_alt = ['A', 'C', 'G']
    num_allels = 0
    alt_read_count = 0
    cum_array = numpy.zeros(12)
    alts = []
    for i in possib_alt:
        row_list = row[i].split(':')
        if int(row_list[1]) > 0:
            num_allels += 1
            allel_count = float(row_list[1])
            alt_read_count += allel_count
            alts += row_list[:2]
            cum_array += allel_count * numpy.array(row_list[2:], dtype=float)
    if alt_read_count > 0:
        cum_array = cum_array / alt_read_count
    clmns_detail = ['base',
                    'count',
                    'avg_mapping_quality',
                    'avg_base_quality',
                    'avg_se_mapping_quality',
                    'num_plus_strand',
                    'num_minus_strand',
                    'avg_pos_as_fraction',
                    'avg_num_mismatches_as_fraction',
                    'avg_sum_mismatch_qualities',
                    'num_q2_containing_reads',
                    'avg_dist_to_q2_start_in_q2_reads',
                    'avg_clipped_length',
                    'avg_dist_to_effective_3p_end']
    res = pandas.Series(ref_split +
                        [num_allels, '_'.join(alts), alt_read_count] +
                        cum_array.tolist() +
                        indel_split,
                        ['REF_'+ x for x in clmns_detail] +
                        ['num_allels'] +
                        ['ALT_'+ x for x in clmns_detail] +
                        ['INDEL_'+ x for x in clmns_detail])
    return res


def addSuffix(x, sfx):
    return [i + sfx for i in x]


def splitVarId(x):
    x_spl = x.split('_')
    if len(x_spl) > 3:
        x_spl = ['_'.join(x_spl[:-2])] + x_spl[-2:]
    return pandas.Series(x_spl, ['ind_id', 'CHROM', 'POS'])


def splitAlleles(x, n_allels=1):
    res = x.split('_')
    # print(res)
    col_names = ['REF', 'ref_DP']
    if len(res) == 3:
        res = [res[0], int(res[1]), '.', 0, int(res[1])]
        return pandas.Series(res, col_names + ['ALT', 'alt_DP', 'DP'])
    for i in range(len(res)/2)[1:]:
        col_names += ['ALT%s' % str(i), 'alt%s_DP' % str(i)]
    DP = sum(map(int, res[1::2]))
    df = pandas.DataFrame({'counts': map(int, res[3::2]),
                           'alleles': res[2::2]})
    df.sort_values(by='counts', ascending=False, inplace=True)
    # print(df)
    alt_DP = df.counts.sum()
    # ALT = ','.join(df.alleles)
    ALT = df.alleles.iloc[0]
    res += [ALT, alt_DP, DP]
    # print(res)
    col_names += ['ALT', 'alt_DP', 'DP']
    # print(col_names)
    return pandas.Series(res, col_names)
#    return pandas.Series(res[:(n_allels + 1) * 2] + res[-1:],
#                         col_names[:(n_allels + 1) * 2] + col_names[-1:])


def mergeClmnsToInfo(df, clmn_list=[]):
    if len(clmn_list) == 0:
        clmn_list = df.columns
    x = None
    for i in clmn_list:
        df['col_i'] = i + '='
        if x is None:
            x = df['col_i'] + df[i].astype(str) + ';'
        else:
            x += df['col_i'] + df[i].astype(str) + ';'
    return x


def writePredAsVcf(pred_df, outp_file, min_DP=0):
    pred_df.reset_index(inplace=True, drop=True)
    # res1 = pred_df.var_id.apply(splitVarId)
    res2 = pred_df.alleles_of.apply(splitAlleles)
    res2 = res2[['ref_DP', 'alt_DP', 'DP']]
    # x = pred_df.merge(res1, left_index=True, right_index=True)
    x = pred_df.merge(res2, left_index=True, right_index=True)
    print 'shape all DP:', pred_df.shape
    x = x[(x['DP_of'] >= min_DP) &
          (x['DP_fa'] >= min_DP) &
          (x['DP_mo'] >= min_DP)]
    print 'shape DP >= ', min_DP
    print x.shape
    required_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    for i in required_fields:
        if i not in x.columns:
            x[i] = '.'
    info_col = [i for i in x.columns if i not in required_fields]
    x['INFO'] = mergeClmnsToInfo(x, info_col)
    x[required_fields + ['INFO']].to_csv(outp_file, sep='\t', index=False,
                                         header=False)
    return 0


def writeTableAsVcf(df, outp_file):
    """write a pandas data frame containing
    CHROM, POS, ID, REF, ALT
    to a vcf file. Remaining columns will preserved in INFO column
    """
    df.reset_index(inplace=True, drop=True)
    required_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    for i in required_fields:
        if i not in df.columns:
            df[i] = '.'
    info_col = [i for i in df.columns if i not in required_fields]
    df['INFO'] = mergeClmnsToInfo(df, info_col)
    df[required_fields + ['INFO']].to_csv(outp_file, sep='\t', index=False,
                                          header=False)
    return 0


def getFieldFromVCF(row, ped_obj, field=6):
#    print row
    ind_id = str(row['ind_id'])
    vcf = ped_obj.getChildsVCF(ind_id)
    cat = 'bcftools view'
#    if os.path.splitext(vcf)[1] == '.gz':
#        cat = 'zcat '
    chrom = str(row['CHROM'])
    pos = str(row['POS'])
#    print ind_id, vcf, cat, chrom, pos
    cmd = ' '.join([cat, vcf, ':'.join([chrom, pos]), '| grep -v ^# | grep ', str(pos)])
    # print cmd
    res = runInShell(cmd, return_output=1)
    if type(res) == int:
        return None
    return res.split('\t')[field]


# def summPred(pred_file, cut_off):
#     tst = train.TrainTest('x',
#                           '/mnt/xfs1/home/asalomatov/projects/variants/variants/ssc_wes_features_noINDEL_noDP.txt',
#                           ['status'],
#                           ['descr'])
#     if not os.path.isfile(pred_file):
#         sys.exit('No such file ; ' + pred_file)
#     df = pandas.read_csv(pred_file)
#     df = df[~df.test_var_id.duplicated()]
#     tst.test_set_y = df['test_labels']
#     # tst.pred_y = df['pred_labels']
#     tst.pred_y = (df['pred_prob'] > cut_off).astype(int)
#     tst.pred_y_prob = df['pred_prob']
#     tst.getMetrics()
#     tst.perf_mertics['method'] = os.path.basename(pred_file)
#     tst.perf_mertics['prob_cutoff'] = cut_off
#     return tst
