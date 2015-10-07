#!/usr/bin/python
'''
look up sequence info for SSCexome.txt  coordinates using pysam, write as headerless vcf.
'''
import pysam, sys

#if len(sys.argv) != 4:
#    print 'usage:'
#    print sys.argv[0], '/pathto.ref.fa /path/to/input.bed /path/to/output.vcf'
#    sys.exit(1)
#
#genref, iosfile, outvcf = sys.argv[1:]
genref = sys.argv[1]

def transformIosLine(x, ref = genref):
    '''
    x is a list with a tab split line from SSCexome.txt file. Transform it to a minimal vcf line
    retaining all of the information in an infor field.
    '''
    fam_id = x[0]
    study = x[1]
    chrom, ios_coord = x[2].split(':')
    if chrom == 'M':
        chrom = 'MT'
    vcf_coord = ios_coord
    variant = x[3]
    #misc_info = x[4:]
    #misc_info.append(study)
    #misc_info = '__'.join(map('_'.join, misc_info))
    #misc_info = '__'.join(['_'.join(x) for x in  misc_info])
    var_desc = variant[0:3]
    ios_ref = ''
    ios_alt = ''
    if var_desc == 'sub':
        x = variant.split('->') 
        ios_ref = x[0][-1]
        ios_alt = x[1][0]
    elif var_desc == 'ins' or var_desc == 'del':
        vcf_coord = str(int(ios_coord) - 1)

    vcf_id = '.'
    qual = '.'
    filt = '.'
    allels = extracAllels(chrom, vcf_coord, variant, genref)
    return '\t'.join([chrom, vcf_coord, vcf_id, allels, qual, filt, 'FAM='+fam_id+';'])
#    info = '.'
#    return '\t'.join([chrom, vcf_coord, vcf_id, allels, qual, filt, info])
    
def extracAllels(chrom, vcf_coord, var_descr, genref):
    '''compute alternative allel from description in bed file.
    must be one of sub(C->T), ins(CCCT), del(5)
    '''
    ref_allel = pysam.faidx(genref, chrom+':'+vcf_coord+'-'+vcf_coord)[1].strip()
    if 'sub' in var_descr:
        yy = var_descr.split('->')
        ref = yy[0][-1]
        if ref == ref_allel:
            return '\t'.join([ref, yy[1][0]])
        else:
            print >> sys.stderr, 'ref allels do not match, exiting...'
            print >> sys.api_version, chrom, vcf_coord, var_descr
            sys.exit(1)
    elif 'ins' in var_descr:
        yy = var_descr.split('(')[1]
        return '\t'.join([ref_allel, ref_allel + yy[:-1]])
    elif 'del' in var_descr:
        yy = var_descr.split('(')[1]
        del_len = int(yy[:-1])
        vcf_coord_end = str(int(vcf_coord) + int(del_len))
        ref = pysam.faidx(genref, chrom+':'+vcf_coord+'-'+vcf_coord_end)[1].strip()
        if ref[0] == ref_allel:
            return '\t'.join([ref, ref_allel])
        else:
            print >> sys.api_version, 'ref allels do not match in del, exiting...'
            print >> sys.api_version, chrom, vcf_coord, var_descr
            sys.exit(1)
    else:
        print >> sys.api_version, 'format not found, exiting...'
        sys.exit(1)

for l in sys.stdin:
    l_list = l.strip().split('\t')
    misc_info = ['_'.join(x.split()) for x in l_list]
    misc_info = '__'.join(misc_info)
    misc_info = misc_info.replace(';','___',100)
    vcfline = transformIosLine(l_list)
    print vcfline+'MISC_INFO='+misc_info

#with open(outvcf, 'a') as fout:
#    with open(iosfile, 'r') as ios:
#        for i, l in enumerate(ios):
#            if i < 2:
#                continue
#            l_list = l.split('\t')
#            misc_info = ['_'.join(x.split()) for x in l_list]
#            #print('__'.join(misc_info))
#            vcfline = transformIosLine(l_list)
#            print vcfline+'MISC_INFO='+'__'.join(misc_info)
#            #if i > 10 : sys.exit(1)
#            fout.write(vcfline+'\n')
