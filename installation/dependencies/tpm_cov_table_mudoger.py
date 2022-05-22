#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
"""A script to calculate TPM values for contigs or genes based on count files
TPM values are defined as in Wagner et al (Theory in Biosciences) 2012. 
      rg x rl x 10^6
TPM = --------------
        flg x T
rg: reads mapped to gene g
rl: read length
flg: feature length 
T: sum of rgxrl/flg for all genes

Coverage calculation
COV =  rg x rl
      ---------
         flg
"""
import sys, pandas as pd, argparse, logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def main(args):
    logging.info("Reading sample info")
    sample_info = pd.read_table(args.sample_info, header=None, index_col=0, names=['avg_read_len'])
    logging.info("Reading gene lengths")
    gene_lengths = pd.read_table(args.gene_lengths, header=None, index_col=0, names=['gene_id','gene_length'])

    df = pd.DataFrame()
    df_cov = pd.DataFrame()

    for fn, sample_name in zip(args.coverage_files, args.sample_names):
        logging.info("Calculating TPM for "+ sample_name)
        ## Read counts per gene for sample
        rg = pd.read_table(fn, index_col=0, header=None, names=['gene_id', 'count'])
        ## Intersect with genes in the gene length file
        rg = rg.loc[list(set(gene_lengths.index).intersection(set(rg.index)))]
        gene_lengths = gene_lengths.loc[list(rg.index)]
        ## Average read length for sample
        rl = sample_info.loc[sample_name,'avg_read_len']
        ## Calculate T for sample
        T = rl * rg['count'].divide(gene_lengths['gene_length']).sum()
        ## Calculate TPM for sample
        tpm = ((1e6*rl)/float(T))*(rg['count'].divide(gene_lengths['gene_length']))
        cov = rl * rg['count'].divide(gene_lengths['gene_length'])
        ## Create dataframe
        TPM = pd.DataFrame(tpm,columns=[sample_name])
        COV = pd.DataFrame(cov,columns=[sample_name])
        ## Concatenate to results
        df = pd.concat([df,TPM],axis=1)
        df_cov = pd.concat([df_cov,COV],axis=1)
    ## Write to file
    df.to_csv(sample_name+".tpm", sep='\t')
    df_cov.to_csv(sample_name+".cov", sep='\t')
    logging.info("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--sample_names', nargs='*', 
            help="Sample names, in the same order as coverage_files")
    parser.add_argument('-c', '--coverage_files', nargs='*', 
            help="Coverage files with tab separated values: 'sequence id, count'")
    parser.add_argument('-i', '--sample_info', 
            help="Tab separated values 'sample_id', 'avg_read_length'")
    parser.add_argument('-l', '--gene_lengths',
            help="Gene lengths in a tsv file")
    args = parser.parse_args()
    main(args)


