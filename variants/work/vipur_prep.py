import pandas, os, sys
from variants import func
import pkg_resources


mis=pandas.read_csv('/mnt/xfs1/home/asalomatov/projects/spark/denovo/SF_denovo_b1-2_b3_b4_b5_mis.csv')
print mis.head()
header_name = os.path.abspath(pkg_resources.resource_filename(
    'variants',
    'header_extra.txt'))
print header_name
dir_name = os.path.dirname(header_name)
print dir_name
vcf_name = '/mnt/xfs1/home/asalomatov/projects/spark/denovo/spark_b1-5_mis.vcf'
cmd = ' '.join(['bash',
                os.path.join(dir_name, 'work', 'mis_for_vipur.sh'),
                vcf_name,
                dir_name])
print cmd
func.runInShell(cmd)


#func.writeTableAsVcf(mis, vcf_name + '.temp')
#func.runInShell('cp ~/projects/variants/variants/header_extra.txt ' +
#                vcf_name)
#func.runInShell('cat ' + 
#                vcf_name + '.temp | sort -V -k1,1 -k2,2 >> ' +
#                vcf_name)
#func.runInShell('rm ' + vcf_name + '.temp')


