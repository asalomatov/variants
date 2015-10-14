import variants
import ped
import func
import pandas as pd
import numpy as np
import sys

class FeaturesVcf:
    """Given list of variants extract features from vcf files"""
    def __init__(self, ped_obj, variants_file ):
        """variants_file contains variants to be annotated with the features from 
        corresponding vcf files, it must be a tab delimeted file containing the following columns
        'ind_id'  'chr' 'pos' 'ref' 'alt', all additional columns will be ignored.
        Example:
        
        ind_id  chr pos ref alt status  descr   vartype
        14075.p1    1   3519049 AC  A   Y   both    ins
        14338.s1    1   6653682 TCC T   Y   Iossifov    ins
        """
        self.ped = ped_obj
        self.variants = pd.read_table(variants_file, index_col=False)

    def _fileHasVariant(self, fname, fam_id, ind_id, chrom, pos_start, pos_end):
        if fname is None:
            return None
        myvars = variants.Variants(fname, fam_id, chrom, pos_start, pos_end)
        try:
            myvars.readFromVcf()
        except:
            return None
        res = self._renameColumns(myvars, ind_id)
        return res


    def _renameColumns(self, myvars, ind_id, num_samples_to_keep=3):
        df = myvars.variants
        if num_samples_to_keep != 1 and num_samples_to_keep != 3:
            return df
        samples_to_keep = {ind_id: 'offspring'}
        if num_samples_to_keep == 3:
            fam_id = self.ped.getFamily(ind_id)
            samples_to_keep[self.ped.getFather(fam_id)] = 'father'
            samples_to_keep[self.ped.getMother(fam_id)] = 'mother'
        a = set(samples_to_keep.keys())
        b = set(myvars.samples)
        if a.intersection(b) != a:
            print 'some amples from ', a, ' are not found in ', b 
            raise
        samples_to_drop = b.difference(a)
        clmns_to_drop = []
        for smpl in samples_to_drop:
            for c in df.columns:
                if smpl in c:
                    clmns_to_drop.append(c)
                else:
                    continue
        if clmns_to_drop:
            df.drop(clmns_to_drop, axis=1, inplace=True)
        #now rename remaining format columns
        new_col_names = list(df.columns)
        for i, c in enumerate(df.columns):
            for key, item in samples_to_keep.iteritems():
                if key in c:
                    new_col_names[i] = c.replace(key, item)
                else:
                    continue
        df.columns = new_col_names
        df['pheno'] = self.ped.isAffected(ind_id)
        df['ind_id'] = ind_id
        return df


    def extractFeatures(self, num_samples_to_keep=3):
        """Given sample id, coordinated, and allels, extract features
        from vcf file
        num_of_samples_to_keep =
        1 - keep only ind_id sample,
        3 - trio, ind_id plus parents
        else - keep all samples without remaning
            """
        res = []
        for i, row in self.variants.iterrows():
            try:
                ind_id = row['ind_id']
                fname = self.ped.getIndivVCF(ind_id)
                if fname is None:
                    continue
                fam_id = self.ped.getFamily(ind_id)
                chrom = row['chr'] 
                pos_vcf = row['pos']
                print ind_id, fam_id, chrom, pos_vcf - 1, fname
                myvars = variants.Variants(fname, fam_id, chrom, pos_vcf - 1, pos_vcf)
                myvars.readFromVcf()
                df = myvars.variants
                if df is None:
                    continue
                if df.empty:
                    continue
                if num_samples_to_keep != 1 and num_samples_to_keep != 3:
                    res.append(df)
                    continue
                samples_to_keep = {ind_id: 'offspring'}
                if num_samples_to_keep == 3:
                    samples_to_keep[self.ped.getFather(fam_id)] = 'father'
                    samples_to_keep[self.ped.getMother(fam_id)] = 'mother'
                a = set(samples_to_keep.keys())
                b = set(myvars.samples)
                if a.intersection(b) != a:
                    print 'some amples from ', a, ' are not found in ', b 
                    raise
                samples_to_drop = b.difference(a)
                clmns_to_drop = []
                for smpl in samples_to_drop:
                    for c in df.columns:
                        if smpl in c:
                            clmns_to_drop.append(c)
                        else:
                            continue
                if clmns_to_drop:
                    df.drop(clmns_to_drop, axis=1, inplace=True)
                #now rename remaining format columns
                new_col_names = list(df.columns)
                for i, c in enumerate(df.columns):
                    for key, item in samples_to_keep.iteritems():
                        if key in c:
                            new_col_names[i] = c.replace(key, item)
                        else:
                            continue
                df.columns = new_col_names
                df['pheno'] = self.ped.isAffected(ind_id)
                df['ind_id'] = ind_id
                res.append(df)
            except:
                continue

      #  x = [i for i in res if not i is None]
      #  y = [i for i in x if not i.empty]
#        df = self.variants.apply(lambda row: self._fileHasVariant(\
#                self.ped.getIndivVCF(row['ind_id']),\
#                self.ped.getFamily(row['ind_id']),\
#                row['chr'], row['pos'] - 1, row['pos']), axis=1)
#        df = df[df.notnull()]
#        c1 =df.apply(lambda x: x.tolist() != [])
#        df = df[c1]
#        res = [x for y in df for x in y.tolist()]
        return res
#
#
#
#        res = []
#        for smpl_id in ['11006.p1']:
#            fam_id = self.ped.getFamily(smpl_id)
#            chrom = '4'
#            pos_start = 115997958
#            pos_end = 115997959
#            if keep_trio and len(self.ped.getParents(fam_id)) < 2:
#                print 'skipping incomplete trio', fam_id  
#                continue
#            fname = self.ped.getIndivVCF(smpl_id)
#            print fname
#            myvars = variants.Variants(fname, fam_id, chrom, pos_start, pos_end)
#            myvars.readFromVcf()
#            if myvars.variants is not None:
#                res.append(myvars.variants)
#        return res
#



if __name__ == '__main__':
    variants = reload(variants)
    ped = reload(ped)
    infile_train = "/mnt/scratch/asalomatov/data/SSCdeNovoCalls/ssc_exome_verified_snp.txt"
    infile_ped = '/mnt/scratch/asalomatov/data/SSCped/SSC.ped'

    myped = ped.Ped(infile_ped, ['collection'])
    myped.getParents(11006)
    myped.addVcf()
    myped.ped.head(30)

    ssc_known = FeaturesVcf(myped, infile_train)
    ssc_known.variants[ssc_known.variants['ind_id'] == '11006.p1']



#    infile_ped = '/mnt/scratch/asalomatov/data/columbia/pcgc_ped.txt'
#    myped = ped.Ped(infile_ped, [])
#    myped.getParents('1-00034')
#    myped.getFather('1-00034')
#    myped.
#    myped.ped.head()
#    myped.ped.shape
#    myped.addVcf(file_pat = '/mnt/scratch/asalomatov/data/columbia/vcf/deco/%s_%s-02_%s-01.annotated-deco.vcf.gz')
#    sum(myped.ped.vcf.notnull())

#    infile_vcf = '/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/rerun200fam/11006-HC-pm50-ann.vcf.gz'
#infile_vcf = '/mnt/scratch/asalomatov/data/columbia/vcf/deco/1-03173_1-03173-02_1-03173-01.annotated-deco.vcf.gz'
#    '/mnt/scratch/asalomatov/data/columbia/vcf/'
#    myvars = variants.Variants(infile_vcf, '11006', myped.ped)
#    myvars.readFromVcf()



#    myvars.samples
#    record = myvars.vcf_reader.next()
#    record.samples
#    myvars._colNamesFormat()
