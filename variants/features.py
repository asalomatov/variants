import variants
import ped
import func
import pandas as pd
import numpy as np
import sys, os
import tempfile
from multiprocessing import Pool

def multi_wrap(args):
    return func.runBamReadcountsRegions(*args)
        
class Features:
    """Given list of variants extract features from vcf files"""
    def __init__(self, ped_obj, variants_file ):
        """variants_file contains variants to be annotated with the features 
        from corresponding vcf files, it must be a tab delimeted file 
        containing the following columns
        'ind_id'  'chr' 'pos' 'ref' 'alt', 'status'. 
        Additional columns may be present.
        ind_id must correspond to samle names in your vcf file(s).
        ped_obj.ped should contain columns with paths to vcf and bam files.
        
        Example variants_file:
        
        ind_id  chr pos ref alt status  descr   vartype
        14075.p1    1   3519049 AC  A   Y   both    ins
        14338.s1    1   6653682 TCC T   Y   Iossifov    ins
        """
        self.ped = ped_obj
        self.family_id = ''
        self.sample_id = ''
        self.father_id = ''
        self.mother_id = ''
        self.sample_vcf = ''
        self.sample_bam = ''
        self.father_bam = ''
        self.mother_bam = ''
        self.is_affected = None
        self.trio_initialized = False
        self.verified_variants = pd.read_table(variants_file, index_col=False)
        self.test_set = pd.DataFrame()
        self.train_set = pd.DataFrame()
        self.sample_features = ''
        self.father_features = ''
        self.mother_features = ''
        self.family_features = ''
    def initTrioFor(self, sample_id):
        """sample_id of a proband or a sibling. Bam files for the trio, and 
        a vcf file for the child have to be specified in ped object.
        """
        self.sample_id = sample_id
        if self.__trioChecked():
            self.trio_initialized = True
            return True
        else:
            self.trio_initialized = False
            return False
        
    def __trioChecked(self):
        """Check if all members, and their files are defined in ped.
        """
        self.family_id = self.ped.getFamily(self.sample_id)       
        if self.family_id is None:
            sys.stderr.write('no family for ' + self.sample_id + '\n')
            return False
        self.father_id = self.ped.getFather(self.family_id)       
        if self.father_id is None:
            sys.stderr.write( 'no father in family ' + self.family_id + '\n')
            return False
        self.mother_id = self.ped.getMother(self.family_id)       
        if self.mother_id is None:
            sys.stderr.write( 'no mother in family ' + self.family_id + '\n')
            return False
        self.sample_vcf = self.ped.getIndivVCF(self.sample_id)       
        if self.sample_vcf is None:
            sys.stderr.write( 'no vcf file for ' + self.sample_id + '\n')
            return False
        self.sample_bam = self.ped.getIndivBAM(self.sample_id)       
        if self.sample_bam is None:
            sys.stderr.write( 'no bam file for ' + self.sample_id + '\n')
            return False
        self.father_bam = self.ped.getIndivBAM(self.father_id)       
        if self.father_bam is None:
            sys.stderr.write( 'no bam file for ' + self.father_id + '\n')
            return False
        self.mother_bam = self.ped.getIndivBAM(self.mother_id)       
        if self.mother_bam is None:
            sys.stderr.write( 'no bam file for ' + self.mother_id + '\n')
            return False
        self.is_affected = self.ped.isAffected(self.sample_id)       
        if self.is_affected is None:
            sys.stderr.write( 'no phenotype for ' + self.sample_id + '\n')
            return False
        return True
    
    def removeTmpDir(self):
        tmpdir = os.path.dirname(self.sample_features)
        sys.stdout.write('removing ' + tmpdir)
        func.runInShell('rm -rf ' + tmpdir)

    def extractFeatures(self):
        """For the defined sample extract variant loci from the vcf file.
        """
        vrs = variants.Variants(self.sample_vcf, self.family_id)
        vrs.readVcfToDF()
        vrs.removeHomRef(self.sample_id)
        vrs.removeHomVar(self.father_id)
        vrs.removeHomVar(self.mother_id)
        temp_dir = tempfile.mkdtemp()
        reg_file = tempfile.mktemp(dir=temp_dir, suffix='_' + self.sample_id \
        +'.region' )
        print reg_file        
        sys.stdout.flush()        
        vrs.vcfDF2regions(reg_file)
        self.sample_features = os.path.join(temp_dir, self.sample_id + '.features')
        print self.sample_features
        self.father_features = os.path.join(temp_dir, self.father_id + '.features')
        print self.father_features
        self.mother_features = os.path.join(temp_dir, self.mother_id + '.features')
        print self.mother_features
        sys.stdout.flush()        

        pool = Pool(3)
        results = pool.map(multi_wrap, \
            [(reg_file, self.sample_bam, self.sample_features), \
            (reg_file, self.father_bam, self.father_features), \
            (reg_file, self.mother_bam, self.mother_features)])
        print results
        if max(results) > 0:
            sys.stderr.write('Feature extraction failed')
            return 1
        return 0
        
        
        
