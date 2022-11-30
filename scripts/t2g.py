#!/usr/bin/env python

import sys, argparse

def create_transcript_list(input, use_name = True, use_version = True):
    r = {}
    for line in input:
        if len(line) == 0 or line[0] == '#':
            continue
        l = line.strip().split(' ')
        if len(l) == 1:
            continue

        if l[1] == 'cdna':
            
            this_transcript = l[0].strip('>')
            this_gene = l[3].split(":")[1].split(".")[0]
            r[this_transcript]=this_gene
    return r



def print_output(output, r, use_name = True):
    f = open(output, "w")
    f.write("transcript_id"+"\t"+"gene_id"+"\n")
    for tx, gene in r.items():
        f.write(tx+"\t"+gene+"\n")
    f.close()


if __name__ == "__main__":


    parser = argparse.ArgumentParser(add_help=True, description='Creates transcript to gene info from GTF files\nreads from standard input and writes to standard output')
    parser.add_argument('--fasta', '-f', dest="fasta")
    parser.add_argument('--output', '-o', dest="output")
    parser.add_argument('--use_version', '-v', action='store_true', help='Use version numbers in transcript and gene ids')
    parser.add_argument('--skip_gene_names', '-s', action='store_true', help='Do not output gene names')
    args = parser.parse_args()



    input = open(args.fasta)

    r = create_transcript_list(input.readlines(), use_name = not args.skip_gene_names, use_version = args.use_version)
    print_output(args.output, r)
