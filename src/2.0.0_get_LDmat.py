import hail as hl
from hail.linalg import BlockMatrix
import pandas as pd
import subprocess
import argparse
import os

def get_UKBB_LDmat_per_locus(ss_name, ld_mat, window_mb, locus):
    print(f"Processing locus: {locus}")
#read in UKBB data
    ht_idx = hl.read_table(f'gs://nygc-bd2disc2-ukbb-ldmatrix/UKBB.{ld_mat}.ldadj.variant.ht').key_by('locus', 'alleles')
    bm = BlockMatrix.read(f'gs://nygc-bd2disc2-ukbb-ldmatrix/UKBB.{ld_mat}.ldadj.bm')
    #read in ss file per locus and create key by locus alleles
    dir = 'gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/' + ss_name + '/'
    snps_file = dir + 'ss/' + ss_name + '_' + window_mb + 'Mb_' + locus + '.txt'
    snps = hl.import_table(snps_file, delimiter = "\t").key_by('chromosome', 'position')
    snps1 = snps.annotate(variant = snps.chromosome + ":" + snps.position + ":" + snps.allele2 + ":" + snps.allele1)
    snps2 = snps.annotate(variant = snps.chromosome + ":" + snps.position + ":" + snps.allele1 + ":" + snps.allele2)
    snps1 = snps1.transmute(**hl.parse_variant(snps1.variant)).key_by('locus', 'alleles')
    snps2 = snps2.transmute(**hl.parse_variant(snps2.variant)).key_by('locus', 'alleles')

    #Join snps to ht index
    ht_idx1 = ht_idx.join(snps1, 'inner')
    ht_idx2 = ht_idx.join(snps2, 'inner')
    ht_idx = ht_idx1.union(ht_idx2)
    ht_idx.export(dir + 'LD/' +  ss_name + '_' + window_mb + 'Mb_' + locus +'_matched.tsv.bgz')
    idx = ht_idx.idx.collect()

    #filter variants
    bm = bm.filter(idx, idx)

    #write new block matrix
    bm_path = dir + 'LD/UKBB_LDmat_' + ss_name + '_' + window_mb + 'Mb_' + locus  
    bm.write(bm_path + '.bm', force_row_major=True)
    BlockMatrix.export(
        bm_path + '.bm',
        bm_path + '.bgz',
        delimiter=' '
    )

def get_UKBB_LDmat(ss_name, ld_mat, window_mb):
  locus_file =  hl.import_table('gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/' + ss_name + '/' + ss_name + '_leadSNPs.tsv', delimiter = "\t")
  locus_file_pd = locus_file.to_pandas()
  for locus in locus_file_pd['locus']:
    get_UKBB_LDmat_per_locus(ss_name, ld_mat, window_mb, locus)
    
    
def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('--ss_name', action="store", dest='ss_name')
  parser.add_argument('--ld_mat', action="store", dest='ld_mat')
  parser.add_argument('--window_mb', action="store", dest='window_mb')
  args = parser.parse_args()
  
  get_UKBB_LDmat(args.ss_name, args.ld_mat, args.window_mb)
  
    
if __name__ == '__main__':
    main()
    
