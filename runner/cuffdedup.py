#!/usr/bin/env python

import sys
import os
import subprocess
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser("Fix cufflinks")
    parser.add_argument('--gtf', dest='gtffile')
    parser.add_argument('--output', dest='outfile')

    return parser

def get_gtf_trans_info(gtf_filename):
    gtf_file = open(gtf_filename)
    trans = {}
    trans_exon=defaultdict(lambda: defaultdict(str))

    for e, line in enumerate(gtf_file):
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        cols = line.split('\t')
        if len(cols) < 8:
            continue

        transid=cols[8].split(" ")[3]
        
        
        if cols[2]=="exon":
            exonnumber=cols[8].split(" ")[5]
            trans_exon[transid][exonnumber] = line
            
        if cols[2]=="transcript":
            if transid not in trans.keys():
                trans[transid] = line
            else:
                continue
            
    return trans, trans_exon
    

def main():
    args = parse_args().parse_args()
    trans, trans_exon = get_gtf_trans_info(args.gtffile)
    
    outfile = open(args.outfile,"w")
    for tr in trans.keys():
        outfile.write(trans[tr]+"\n")
        
        for exs in trans_exon[tr]:
            outfile.write(trans_exon[tr][exs]+"\n")

    outfile.close()
    
    
if __name__=="__main__":
    sys.exit(main())
### read in gtf
## line 

## transcript start dict

## transcript end dict