"""Class and container for pedigree information, vcf, and bam file by sample"""

import pandas as pd
import re 
import func

class Ped:
    """Family_ID - '.' or '0' for unknown
    Individual_ID - '.' or '0' for unknown
    Paternal_ID - '.' or '0' for unknown
    Maternal_ID - '.' or '0' for unknown
    Sex - '1'=male; '2'=female; ['other', '0', '.']=unknown
    Phenotype - '1'=unaffected, '2'=affected, ['-9', '0', '.']= missing"""

    def __init__(self, ped_file_name, extra_column_names=[]):
        """read ped file into pandas data frame"""
        self.fname = ped_file_name
        self.ped = pd.read_table(self.fname)
        self.ped.columns = ['fam_id', 'ind_id', 'fa_id', 'mo_id', 'sex', 'pheno'] + extra_column_names
        self.ped.replace(['.', '0', 0, -9, '-9'], [None]*5, inplace=True)

    def addVcf(self, field='fam_id', file_pat='/mnt/ceph/asalomatov/SSC_Eichler/rerun/ssc%s/%s-JHC-vars.vcf.gz'):
        num_subst = len(re.findall('\%s', file_pat))
        print num_subst, ' substitutions found'
        if num_subst > 0:
            x = self.ped[field].apply(lambda f: func.checkFile(file_pat % ((f,) * num_subst)))
            self.ped['vcf'] = pd.Series(x, index=self.ped.index)
        else:
            self.ped['vcf'] = file_pat

    def addBam(self, field='ind_id', file_pat='/mnt/ceph/asalomatov/SSC_Eichler/data_S3/%s-*.bam', num_subst=2):
        pass

    def getAllMembers(self, family_id):
        return self.ped['ind_id'][self.ped['fam_id'] == family_id].tolist()

    def getProbands(self, family_id):
        return self.ped['ind_id'][(self.ped['fam_id'] == family_id) & (self.ped['pheno'] == 2)].tolist()

    def getSiblings(self, family_id):
        return self.ped['ind_id'][(self.ped['fam_id'] == family_id) & (self.ped['pheno'] == 1) \
                & ~self.ped['fa_id'].isnull() & ~self.ped['mo_id'].isnull() ].tolist()

    def getParents(self, family_id):
        return self.ped['ind_id'][(self.ped['fam_id'] == family_id) &  \
                self.ped['fa_id'].isnull() & self.ped['mo_id'].isnull() ].tolist()

    def getFather(self, family_id):
        res = self.ped['ind_id'][(self.ped['fam_id'] == family_id) & (self.ped['sex'] == 1) &  \
                self.ped['fa_id'].isnull() & self.ped['mo_id'].isnull() ]
        assert len(res) == 1
        return res.iloc[0]

    def getMother(self, family_id):
        res = self.ped['ind_id'][(self.ped['fam_id'] == family_id) & (self.ped['sex'] == 2) &  \
                self.ped['fa_id'].isnull() & self.ped['mo_id'].isnull() ]
        assert len(res) == 1
        return res.iloc[0]

    def getChildsFather(self, individial_id):
        res = self.ped['fa_id'][(self.ped['ind_id'] == individial_id)]
        assert len(res) == 1
        return res.iloc[0]

    def getChildsMother(self, individial_id):
        res = self.ped['mo_id'][(self.ped['ind_id'] == individial_id)]
        assert len(res) == 1
        return res.iloc[0]

    def isAffected(self, individial_id):
        res = self.ped['pheno'][(self.ped['ind_id'] == individial_id)] == 2
        assert len(res) == 1
        return res.iloc[0]

    def getIndivVCF(self, individial_id):
        res = self.ped['vcf'][(self.ped['ind_id'] == individial_id)]
        assert len(res) == 1
        return res.iloc[0]

    def getIndivBAM(self, individial_id):
        res = self.ped['bam'][(self.ped['ind_id'] == individial_id)]
        assert len(res) == 1
        return res.iloc[0]

    def getFamily(self, individial_id):
        res = self.ped['fam_id'][(self.ped['ind_id'] == individial_id)]
        assert len(res) == 1
        return res.iloc[0]

    def getFamilyVCF(self, family_id):
        res = self.ped['vcf'][(self.ped['fam_id'] == family_id)]
        res = res.unique()
        assert len(res) == 1
        return res[0]

    def getFamilyBam(self, family_id):
        res = self.ped['bam'][(self.ped['fam_id'] == family_id)]
        res = res.unique()
        assert len(res) == 1
        return res[0]

    def getAllProbands(self):
        res = self.ped['ind_id'][self.ped['pheno'] == 2]
        res = res.tolist()
        return res

    def getAllTrios(self):
        fam = self.ped['fam_id'].unique()
        res = [x for x in fam if len(self.getAllMembers(x)) == 3]
        return res

    def getAllQuads(self):
        fam = self.ped['fam_id'].unique()
        res = [x for x in fam if len(self.getAllMembers(x)) == 4]
        return res

    def isTrio(self, family_id):
        res = len(self.ped['fam_id'][(self.ped['fam_id'] == family_id)]) == 3
        return res

    def isQuad(self, family_id):
        res = len(self.ped['fam_id'][(self.ped['fam_id'] == family_id)]) == 4
        return res

if __name__ == '__main__':
    infile = '/mnt/scratch/asalomatov/data/SSCped/SSC.ped'
    myped=Ped(infile, ['collection'])
    myped.addVcfSSC()
