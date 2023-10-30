from math import remainder
import numpy as np
import pandas as pd
import argparse
import os
import sys
import pickle
from gtfparse import read_gtf

def parse_args():
    parser = argparse.ArgumentParser("Read kallisto")
    parser.add_argument('--dir', dest='dir')
    parser.add_argument('--outputfile', dest='outputfile')
    parser.add_argument('--pattern', dest='pattern')
    parser.add_argument('--type', dest='type', default='TPM')
    parser.add_argument('--gtf', dest='gtf')
    parser.add_argument('--meta', dest='meta')
    return parser

if __name__=="__main__":
    args = parse_args().parse_args()
    print(args.pattern)
    reads = []
    for x in os.listdir(args.dir):
        print(x)
        #if args.pattern in x:
        tmp = pd.read_csv(os.path.join(args.dir, x, "abundance.tsv"), sep="\t")
        print(tmp.head())
        if args.type=='TPM':
            reads.append(pd.Series(tmp['tpm'], name=x))
        else:
            reads.append(pd.Series(tmp['est_counts'], name=x))
            #reads.append(pd.Series(tmp['TPM'],name=x)) 
    #reads.append(pd.Series(tmp['Name'], name="feature_id"))
    finaldf = pd.concat(reads, axis=1)
    print(finaldf.head())
    print(finaldf.shape)
    
    ## read generate tx2gene
    gtf = pd.read_csv(args.gtf, sep='\t')
    print(gtf.head())
    ###check salmon out put id version
  
    genemapper = dict(zip(gtf['transcript_id'], gtf['gene_id']))
    #pickle.dump(genemapper, open("./gene_mapper.pkl", "wb"))
    
    #genemapper = pickle.load(open("./gene_mapper.pkl", "rb"))
    #
    
    ### make a file for iso_ktsp 

    
    meta = pd.read_csv(args.meta, sep="\t")
    meta['sample_id'] = list(map(lambda x: x.replace(".sra", ""), meta['sample_id']))

    print(finaldf.loc[:, np.isin(finaldf.columns, meta['sample_id'])])
    #keep only columns that are in metadata
    finaldf = finaldf.loc[:, np.isin(finaldf.columns, meta['sample_id'])]

    # map_group = {x:tmp[e] for e,x in enumerate(meta['group'].unique())}
    meta['name_isoktsp'] = meta['sample_id']+"_"+meta['group'].astype('str')
    print(meta.head())
    isodf = finaldf.copy()
    isodf.columns = meta['name_isoktsp'].to_list() 
    isodf['isoform_id'] = tmp['target_id']
    isodf = isodf[['isoform_id']+meta['name_isoktsp'].to_list()]

    ##isoform_id

    def mapping(x):
        try:
            return genemapper[x]
        except:
            return "?"
    gene_id = list(map(lambda x: mapping(x), tmp['target_id']))

    isodf['isoform_id'] = pd.Series(gene_id) + "," + isodf['isoform_id']
    isodf.columns = [""] + meta['name_isoktsp'].to_list() 

    isodf.to_csv(os.path.join(args.outputfile, "raw_kal_isoktsp_count.txt"), sep="\t", index=False)
    finaldf['gene_id'] = gene_id
    finaldf['feature_id'] = tmp['target_id']
    finaldf.to_csv(os.path.join(args.outputfile, "kal_count.csv"), index=False, sep="\t")

    f=open(os.path.join(args.outputfile, "raw_kal_isoktsp_count.txt"), "r")
    firstline, rest = f.readline(), f.read()
    t=open(os.path.join(args.outputfile, "kal_isoktsp_count.txt"), "w")
    t.write(firstline[1:])
    t.write(rest)
    t.close()
    f.close()
    print("Done")
