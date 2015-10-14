import pandas as pd
import vcf
import sys, os, re
import func
                                                
class Variants:
    """Class describing a vcf file. Consists of pandas data frame, and metadata found in the vcf header.
       Start, end coordinates are zero-based, half-open """
    def __init__(self, fname, family_id, chrom=None, start=None, end=None):
        self.fname = fname
        self.family_id = family_id
        self.vcf_reader = vcf.Reader(open(self.fname, 'r'))
        if not chrom is None:
            self.vcf_reader = self.vcf_reader.fetch(chrom, start, end)
        self.current_record = None
        self.variants = pd.DataFrame()
        self.chrom = chrom
        self.start = start
        self.end = end
        self.contigs = self.vcf_reader.contigs
        self.filters = self.vcf_reader.filters
        self.formats = self.vcf_reader.formats
        self.infos = self.vcf_reader.infos
        self.metadata = self.vcf_reader.metadata
        self.samples = self.vcf_reader.samples
        self.samples_to_keep = self.vcf_reader.samples
        self.required_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
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
            print k, f.type, f.desc

    def describeFormatFields(self):
        for k, f in self.formats.items():
            print k, f.type, f.desc

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
            #print 'processing sample ', smpl
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
           # print '  '
           # print 'parsing sample ', smpl.sample
            if not smpl.called:
                line_format += [None] * (len(self.format_fields) + 3 +1)
                continue
            else:
                for format_f in self.format_fields:
           #         print line_format
           #         print 'parsing', format_f
                    try:
                        format_v = getattr(smpl.data, format_f, None)
                        if format_f == 'AD':
                            if format_v is None:
                                line_format += [None] * 2
                            else:
                                line_format += format_v
                                #print format_f, format_v, line_format
                        elif format_f == 'PL':
                            if format_v is None:
                                line_format += [None] * 3
                            else:
                                line_format += format_v
                                #print format_f, format_v, line_format
                        else:
                            line_format.append(format_v)
                            #print format_f, format_v, line_format
                    except:
                        print self.current_record, self.current_record.INFO, smpl, smpl.called
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
            #print self.current_record, self.current_record.INFO, self.current_record.samples
            #print ' '
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
               # print line_info
                #FORMAT
                line_format = self._parseFormat(i_alt)
               # print line_format
                line_extra = [self._varType(), record.is_transition]
               # print line_extra
                lines.append(line_required + line_info + line_format + line_extra)
        #        print line
        #        if i > 2:
        #            sys.exit()
        if not lines:
            print 'no variants in this region'
        else:
            self.variants = pd.DataFrame(lines, dtype=str)
            self.variants.columns = self.required_fields + col_names_info  + col_names_format + col_names_extra
            self.variants['family_id'] = self.family_id
            #self.variants[1:] = self.variants[1:].convert_objects(convert_numeric=True)
            print "%s multiallelic sites skipped" % num_multiallelic
    
    def keepOnlyPossibleDonovos():
        pass

    def saveAsVcf(self):
        pass

if __name__ == '__main__':
    pass
