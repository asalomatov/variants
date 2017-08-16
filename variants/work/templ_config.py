from variants import func


'''
alt_allele_frac_range:
- 0.2
- 0.8
bam_pattern: /mnt/xfs1/scratch/asalomatov/data/SPARK/bam/%s.bam
bam_readcount: /mnt/xfs1/scratch/asalomatov/software/installs/bin/bam-readcount
db_nsfp:
  M_CAP: 0.025
  M_CAP_pred:
  - D
  cadd_phred: 25
  combined:
    cadd_phred: 15
    polyphen2_pred:
    - D
    sift_pred:
    - D
  metaSVM_pred:
  - D
denovo_definition: 1
genome_build: 37
genome_ref: /mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
known_variants: null
max_cohort_freq: 1
min_DP: 10
output_directory: /mnt/ceph/users/asalomatov/spark/denovo/%(dr)s/snp/%(clr)s/%(b)s
ped_file: null
ped_file_extended: /mnt/xfs1/scratch/asalomatov/data/SPARK/ped/spark_spID_%(b)s_ext_%(clr)s.ped
population_AF: 0.01
snpeff:
  biotype:
  - protein_coding
  effect_dmgmis:
  - missense_variant
  effect_lof:
  - exon_loss_variant
  - frameshift_variant
  - stop_gained
  - stop_lost
  - start_lost
  - splice_acceptor_variant
  - splice_donor_variant
  - splice_region_variant
  effect_synon:
  - synonymous_variant
  genes:
  - ADNP
  - ADSL
  - AHDC1
  - ALDH5A1
  - ANK2
  - ANKRD11
  - ARID1B
  - ARX
  - ASH1L
  - ASXL3
  - AUTS2
  - BCKDK
  - BCL11A
  - CACNA1C
  - CDKL5
  - CHD2
  - CHD7
  - CHD8
  - CREBBP
  - DDX3X
  - DHCR7
  - DMPK
  - DSCAM
  - DYRK1A
  - EHMT1
  - EP300
  - FMR1
  - FOXP1
  - GIGYF2
  - GRIN2B
  - IQSEC2
  - KIAA2022
  - MBD5
  - MBOAT7
  - MECP2
  - MED13L
  - MYT1L
  - NCKAP1
  - NF1
  - NHE6
  - NIPBL
  - NLGN2
  - NLGN3
  - NRXN1
  - NRXN2
  - NRXN3
  - NSD1
  - PACS1
  - POGZ
  - POMGNT1
  - PTCHD1
  - PTEN
  - RAI1
  - RIMS1
  - SCN1A
  - SCN2A
  - SETBP1
  - SETD2
  - SETD5
  - SHANK2
  - SHANK3
  - SLC6A1
  - SLC9A6
  - SON
  - STXBP1
  - SUV420H1
  - SYNGAP1
  - TBCK
  - TRIO
  - TRIP12
  - TSC1
  - TSC2
  - TSHZ3
  - UBE3A
  - UPF3B
  - VPS13B
  - WAC
  - ZBTB20
  impact_dmgmis: []
  impact_lof:
  - HIGH
target_bed: /mnt/xfs1/scratch/asalomatov/data/SPARK/info/HG19_vcrome2.1_with_PKv2_pm50.bed
variant_type: SNP
vcf_pattern: null
'''

t = func.readYml('/mnt/xfs1/home/asalomatov/projects/variants/variants/config_files/cfg.yml')
denovo_def = 3
outp_dir = 'def' + str(denovo_def)

for vartype in ['snp', 'indel']:
    for clr in ['hc', 'fb', 'pl']:
        for bbbatch in ['b1-2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9']:
            config_file = '/mnt/ceph/users/asalomatov/spark/denovo/%(outp_dir)s/cfg_spark_%(bbbatch)s_%(vartype)s_%(clr)s.yml' % locals()
            print config_file
            t['output_directory'] = '/mnt/ceph/users/asalomatov/spark/denovo/%(outp_dir)s/%(vartype)s/%(clr)s/%(bbbatch)s' % locals()
            t['ped_file_extended'] = '/mnt/xfs1/scratch/asalomatov/data/SPARK/ped/spark_spID_%(bbbatch)s_ext_%(clr)s.ped' % locals()
            t['denovo_definition'] = denovo_def
            t['variant_type'] = vartype.upper()
            print t['output_directory']
            print t['ped_file_extended']
            print t['denovo_definition']
            print t['variant_type']
            func.dumpYml(config_file, t)

#ped_file: null
#ped_file_extended: /mnt/xfs1/scratch/asalomatov/data/SPARK/ped/spark_spID_%(b)s_ext_%(clr)s.ped
