from __future__ import print_function
import pandas as pd
import vcf
import sys, os, re
import func
import collections

class Variants:
    """Class describing a vcf file. Consists of pandas data frame, and metadata found in the vcf header.
       Start, end coordinates are zero-based, half-open """
    def __init__(self, fname, family_id, chrom=None, start=None, end=None):
        self.fname = fname
        self.family_id = family_id
        self.vcf_reader = None
        if not chrom is None:
            self.vcf_reader = self.vcf_reader.fetch(chrom, start, end)
        self.current_record = None
        self.variants = pd.DataFrame()
        self.chrom = chrom
        self.start = start
        self.end = end
        self.required_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']

    def initReader(self):
        self.vcf_reader = vcf.Reader(open(self.fname, 'r'))
        self.contigs = self.vcf_reader.contigs
        self.filters = self.vcf_reader.filters
        self.formats = self.vcf_reader.formats
        self.infos = self.vcf_reader.infos
        self.metadata = self.vcf_reader.metadata
        self.samples = self.vcf_reader.samples
        self.samples_to_keep = self.vcf_reader.samples
        self.info_fields = self.infos.keys()
        self.format_fields = self.formats.keys()
        self.reference = self.metadata['reference']

    def guessCaller(self):
        y = [x.lower() for x in self.metadata.keys()]
        y = ' '.join(y)
        if 'haplotypecaller' in y:
            return 'HaplotypeCaller'
        elif 'freebayes' in y:
            return 'Freebayes'
        elif 'platypus' in y:
            return 'Platypus'
        else:
            return 'Unknown'

    def infoFieldDescr(self, field_name):
        return self.infos[field_name].desc
    
    def infoFieldType(self, field_name):
        return self.infos[field_name].type

    def formatFieldDescr(self, field_name):
        return self.formats[field_name].desc

    def formatFieldType(self, field_name):
        return self.formats[field_name].type

    def describeInfoFields(self):
        for k, f in self.infos.items():
            print(' '.join([k, f.type, f.desc]))

    def describeFormatFields(self):
        for k, f in self.formats.items():
            print(' '.join([k, f.type, f.desc]))

    def _formatFilter(self, filt):
        if not filt:
            return 'PASS'
        else:
            return ','.join(filt)

    def _parseInfo(self, i_alt_allel):
        line_info = []
        for info_f in self.info_fields:
            if not info_f in self.current_record.INFO:
                line_info.append(None)
                continue
            info_v = self.current_record.INFO[info_f]
            if not isinstance(info_v, list):
                line_info.append(info_v)
            elif not info_v:
                line_info.append(None)
            else:
                line_info.append(info_v[i_alt_allel])
        return line_info

    def _sample2member(self, sample_name):
        return sample_name
#        m = re.search('\.(s\d)|\.(p\d)|\.(fa)|\.(mo)', sample_name)
#        if not m:
#            return sample_name
#        else:
#            ml = [i for i in m.groups() if i is not None]
#            if len(ml) > 1:
#                print 'More than one match in sample name'
#                raise
#            else:
#                return ml[0]
#        for m in ['p1', 'p2', 'p3', 'p4', 'p5']:
#            if m in sample_name:
#                return m
#            elif 'mo' in sample_name:
#                return 'mo'
#            else:
#                return 'fa'
    
    def _colNamesFormat(self): 
        col_names_format = []
#        memb = [self._sample2member(x.sample) for x in self.current_record.samples]
        for smpl in self.samples:
            for format_f in self.format_fields:
                    if format_f == 'AD':
                        col_names_format.append('_'.join(['format', smpl, 'ref', format_f])) 
                        col_names_format.append('_'.join(['format', smpl, 'alt', format_f])) 
                    elif format_f == 'PL':
                        col_names_format.append('_'.join(['format', smpl, '0', format_f])) 
                        col_names_format.append('_'.join(['format', smpl, '1', format_f])) 
                        col_names_format.append('_'.join(['format', smpl, '2', format_f])) 
                    else:
                        col_names_format.append('_'.join(['format', smpl, format_f])) 
            col_names_format.append('_'.join([smpl, 'gt_type']))
        return col_names_format

    def _parseFormat(self, i_alt_allel):
        line_format = []
        for smpl in self.current_record.samples:
            if not smpl.called:
                line_format += [None] * (len(self.format_fields) + 3 +1)
                continue
            else:
                for format_f in self.format_fields:
                    try:
                        format_v = getattr(smpl.data, format_f, None)
                        if format_f == 'AD':
                            if format_v is None:
                                line_format += [None] * 2
                            else:
                                line_format += format_v
                        elif format_f == 'PL':
                            if format_v is None:
                                line_format += [None] * 3
                            else:
                                line_format += format_v
                        else:
                            line_format.append(format_v)
                    except:
                        print(' '.join([self.current_record,
                                        self.current_record.INFO,
                                        smpl, smpl.called]))
                        raise
            line_format.append(smpl.gt_type)
        return line_format
    def _varType(self):
        if self.current_record.is_deletion:
            return 'del'
        elif self.current_record.is_indel and not self.current_record.is_deletion:
            return 'ins'
        elif self.current_record.is_snp:
            return 'snp'
        else:
            return self.current_record.var_type
    

    def readFromVcf(self, min_depth=8):
        lines = []
        num_multiallelic = 0
        col_names_info = ['_'.join(['info', x]) for x in self.info_fields]
        col_names_format = self._colNamesFormat()
        col_names_extra = ['vartype', 'is_transition']
        for i, record in enumerate(self.vcf_reader):
            if not record.INFO['DP'] or int(record.INFO['DP']) < min_depth:
                continue
            self.current_record = record
            if len(record.ALT) > 1: #for now skip all multiallelic
                num_multiallelic += 1
                continue
            for i_alt, nucl_alt in enumerate(record.ALT):
                line_required = [record.CHROM, 
                        record.POS, 
                        record.ID, 
                        record.REF, 
                        nucl_alt.sequence, #record.ALT[0].sequence,
                        record.QUAL, 
                        self._formatFilter(record.FILTER)]
                #INFO
                line_info = self._parseInfo(i_alt)
                #FORMAT
                line_format = self._parseFormat(i_alt)
                line_extra = [self._varType(), record.is_transition]
                lines.append(line_required + line_info + line_format + line_extra)
        if not lines:
            print('no variants in this region')
        else:
            self.variants = pd.DataFrame(lines, dtype=str)
            self.variants.columns = self.required_fields + col_names_info  + col_names_format + col_names_extra
            self.variants['family_id'] = self.family_id
            #self.variants[1:] = self.variants[1:].convert_objects(convert_numeric=True)
            print("%s multiallelic sites skipped" % num_multiallelic)
    
    def readVcfToDF(self, chunk_size=None):
        """read vcf file into pandas DF without parsing"""
        # after github/hammerlab/varcode/vcf but keeping sample information
        path = self.fname        
        compression = None
        if path.endswith(".gz"):
            compression = "gzip"
        elif path.endswith(".bz2"):
            compression = "bz2"
        cat = 'cat'
        if compression is not None:
            cat = 'zcat'
        cmd = ' '.join([cat, path, '| grep ^# | grep -v ^##'])
        vcf_clmns = func.runInShell(cmd, True).split('\t')
        vcf_clmns = [x.strip() for x in vcf_clmns]        
        vcf_clmns = [x.strip('#') for x in vcf_clmns]        
        vcf_field_types = collections.OrderedDict()
        for i in vcf_clmns:
            vcf_field_types[i] = str
        vcf_field_types['POS'] = int
        reader = pd.read_table(
            path,
            compression=compression,
            comment="#",
            chunksize=chunk_size,
            dtype=vcf_field_types,
            names=list(vcf_field_types),
            usecols=range(len(vcf_field_types)))
        self.variants = reader
        return 0

    def vcfDF2regions(self, reg_file_name, vartype='SNP'):
        """ Create regions file to supply to bam-readcount.
        CHROM\tSTART\tEND, 1-based
        """
        if vartype.lower() == 'snp':
            ref_len = self.variants.REF.apply(len)
            alt_len_min = self.variants.ALT.apply(lambda x:
                                                  min(map(len, x.split(','))))
            alt_len_max = self.variants.ALT.apply(lambda x:
                                                  max(map(len, x.split(','))))
            # following conditions to capture SNPs and MNPs
            c1 = ref_len == alt_len_min
            c2 = ref_len == alt_len_max
            print('number of MNPs = %s' % sum((ref_len > 1) & (c1 | c2)))
            print('number of SNPs = %s' % sum((ref_len == 1) & (c1 | c2)))
            self.variants['pos_end'] = self.variants['POS'] + ref_len - 1
            self.variants[['CHROM', 'POS', 'pos_end']][c1 | c2].to_csv(
                reg_file_name, sep="\t", header=False, index=False)
        elif vartype.lower() == 'indel':
            ref_len = self.variants.REF.apply(len)
            alt_len_min = self.variants.ALT.apply(lambda x:
                                                  min(map(len, x.split(','))))
            alt_len_max = self.variants.ALT.apply(lambda x:
                                                  max(map(len, x.split(','))))
            # following conditions to capture INDELS
            c1 = (ref_len > 1) & (alt_len_min == 1)  # deletion
            c2 = (ref_len == 1) & (alt_len_max > 1)  # insertion
            print('number of DELETIONS = %s' % sum(c1))
            print('number of INSERTIONS = %s' % sum(c2))
            # deletions should have POS + 1
            self.variants.ix[c1, 'POS'] = self.variants.POS[c1] + 1
            self.variants.ix[:, 'pos_end'] = self.variants.POS
            self.variants[['CHROM', 'POS', 'pos_end']][c1 | c2].to_csv(
                reg_file_name, sep="\t", header=False, index=False)
        else:
            sys.exit('vcfDF2regions: Only SNPs or INDELs are supported')

    def removeHomRef(self, sample_name):
        """Remove non variant loci. Remove loci with missing genotype from
        self.variants DF
        """
        gt = self.variants[sample_name].apply(lambda i: i.split(':')[0])
        gt = [i.strip() for i in gt]
        self.variants[sample_name+'_gt'] = gt        
        c1 = self.variants[sample_name+'_gt'].isin(['0/0'])
        c2 = self.variants[sample_name+'_gt'].apply(lambda i: '.' in i)        
        self.variants = self.variants[(~c1) & (~c2)]         
        return 0

    def removeUndefined(self, sample_name):
        """Remove genotypes like ./. from self.variants DF
        """
        gt = self.variants[sample_name].apply(lambda i: i.split(':')[0])
        gt = [i.strip() for i in gt]
        self.variants[sample_name+'_gt'] = gt        
        c1 = self.variants[sample_name+'_gt'].apply(lambda i: '.' in i)        
        self.variants = self.variants[~c1]         
        return 0

    def removeNaN(self, sample_name):
        """Remove NaN from self.variants DF
        """
        self.variants = self.variants[~self.variants[sample_name].isnull()]
        return 0

    def removeNoGT(self, sample_name):
        """Keep only ['0/0', '0/1', '1/1'] from self.variants DF
        """
        c1 = self.variants[sample_name].str.contains('\d/\d')
        c2 = self.variants[sample_name].str.contains('\d|\d')
        self.variants = self.variants[c1 | c2]
        return 0

    def removeHomVar(self, sample_name):
        """Remove loci where sample_name has 1/1 genotype, or undefined 
        genotype"""
        gt = self.variants[sample_name].apply(lambda i: i.split(':')[0])
        gt = [i.strip() for i in gt]
        self.variants[sample_name+'_gt'] = gt        
        c1 = self.variants[sample_name+'_gt'].isin(['1/1'])
        c2 = self.variants[sample_name+'_gt'].apply(lambda i: '.' in i)        
        c3 = self.variants[sample_name+'_gt'].apply(lambda i: '0' in i)        
        self.variants = self.variants[(~c1) & (~c2) & c3]         
        return 0

    def getFieldFromVCF(row, ped_obj, field=6):
        ind_id = row['ind_id']
        vcf = ped_obj.getIndivVCF(ind_id)
        cat = 'bcftools view'
    #    if os.path.splitext(vcf)[1] == '.gz':
    #        cat = 'zcat '
        chrom = str(row['CHROM'])
        pos = str(row['POS'])
        cmd = ' '.join([cat, vcf, ':'.join([chrom, pos]), '| grep -v ^# | grep ', str(pos)])
        res = func.runInShell(cmd, return_output=1)
        if type(res) == int:
            return None
        return res.split('\t')[field]

    def fileHasVariant(fname, fam_id, chrom, pos_start, pos_end):
        myvars = Variants(fname, fam_id, chrom, pos_start, pos_end)
        next_var = None
        fam_var = []
        for next_var in myvars.vcf_reader:
            alt_allel = []
            for nucl_alt in next_var.ALT: 
                alt_allel.append(nucl_alt.sequence)
            for v in alt_allel:
                fam_var.append('_'.join([fam_id, next_var.CHROM, str(next_var.POS), next_var.REF, v]))
        return fam_var

    def keepOnlyPossibleDenovos(self, smpl_ch, smpl_fa, smpl_mo):
        """child != 0/0 and both parents == 0/0"""
        gt_ch = self.variants[smpl_ch].apply(lambda i: i.split(':')[0].strip())
        #gt = [i.strip() for i in gt]
        self.variants[smpl_ch + '_gt'] = gt_ch        
        gt_fa = self.variants[smpl_fa].apply(lambda i: i.split(':')[0].strip())
        #gt = [i.strip() for i in gt]
        self.variants[smpl_fa + '_gt'] = gt_fa        
        gt_mo = self.variants[smpl_mo].apply(lambda i: i.split(':')[0].strip())
        #gt = [i.strip() for i in gt]
        self.variants[smpl_mo + '_gt'] = gt_mo        
        c1 = self.variants[smpl_ch + '_gt'].isin(['1/1', '0/1', '0/2', '1|1', '0|1', '0|2'])
        c2 = self.variants[smpl_fa + '_gt'].isin(['0/0', '0|0'])
        c3 = self.variants[smpl_mo + '_gt'].isin(['0/0', '0|0'])
        self.variants = self.variants[c1 & c2 & c3]
        return 0

    def saveAsVcf(self):
        pass

if __name__ == '__main__':
    pass
