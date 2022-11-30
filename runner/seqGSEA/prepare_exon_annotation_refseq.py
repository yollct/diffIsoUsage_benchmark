#!/usr/bin/python

####
#### script from https://github.com/Wanvdphelys/SeqGSEA/blob/master/inst/extscripts/prepare_exon_annotation_refseq.py
####
# to prepare exon annotations from RefSeq gene annotation gtf file
# modified from prepare_exon_annotation.py released within DEXSeq R package
# 
# Differences are: (1) group gene fimily together 
# (2) remove overlapping exons of genes not sharing the names
# (3) randomly pick up a location for a muti-located gene (incuding strand) 

import sys, collections, itertools, os.path, operator 
GENE_FAMILY_PREFIX_MIN = 2

def read_through_gene_rescue(gene_list): 
  """
  Return 1 if gene_list contains like: [GENE1-GENE2, GENE2, GENE2]
  Return 0 otherwise
  """
  g_list = []
  min_len = 100
  for g in gene_list : 
    sub_g = g[g.find('-')+1:]
    if min_len > len(sub_g) : 
      min_len = len(sub_g)
    g_list.append(sub_g)
  if len(os.path.commonprefix(g_list)) >= min(GENE_FAMILY_PREFIX_MIN, min_len) :
    return 1
  else :
    return 0

if len( sys.argv ) != 3:
   sys.stderr.write( "Script to prepare annotation for SeqGSEA.\n\n" )
   sys.stderr.write( "Usage: python %s <in.gtf> <out.gff>\n\n" % os.path.basename(sys.argv[0]) )
   sys.stderr.write( "This script takes an annotation file in GTF format\n" )
   sys.stderr.write( "and outputs a 'flattened' annotation file suitable for use\n" )
   sys.stderr.write( "with the count_in_exons.py script.\n" )
   sys.exit(1)

try:
   import HTSeq
except ImportError:
   sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )   
   sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )   
   sys.exit(1)

gtf_file = sys.argv[1]
out_file = sys.argv[2]

# Step 1: Store all exons with their gene and transcript ID 
# in a GenomicArrayOfSets

exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
for f in HTSeq.GFF_Reader( gtf_file ):
   if f.type != "exon":
      continue
   f.attr['gene_id'] = f.attr['gene_id'].replace( ":", "_" )
   exons[f.iv] += ( f.attr['gene_id'], f.attr['transcript_id'] )


# Step 2: Form sets of overlapping genes

# We produce the dict 'gene_sets', whose values are sets of gene IDs. Each set
# contains IDs of genes that overlap, i.e., share bases (on the same strand).
# The keys of 'gene_sets' are the IDs of all genes, and each key refers to
# the set that contains the gene.
# Each gene set forms an 'aggregate gene'.

exon_omit_list = list()
gene_sets = collections.defaultdict( lambda: set() )
for iv, s in exons.steps():
   if len(s) == 0:
     exon_omit_list.append(1)
     continue
   gene_list = list()
   min_len = 100
   for gene_id, transcript_id in s: 
     gene_list.append(gene_id)
     if min_len > len(gene_id):
       min_len = len(gene_id)
   commonprefix_gene_list = os.path.commonprefix(gene_list)
   if len(commonprefix_gene_list) >= min(GENE_FAMILY_PREFIX_MIN, min_len) or read_through_gene_rescue(gene_list) :
     exon_omit_list.append(0)
     # For each step, make a set, 'full_set' of all the gene IDs occuring
     # in the present step, and also add all those gene IDs, whch have been
     # seen earlier to co-occur with each of the currently present gene IDs.
     full_set = set()
     for gene_id, transcript_id in s:
       full_set.add( gene_id )
       full_set |= gene_sets[ gene_id ]
     # Make sure that all genes that are now in full_set get associated
     # with full_set, i.e., get to know about their new partners
     for gene_id in full_set:
       assert gene_sets[ gene_id ] <= full_set
       gene_sets[ gene_id ] = full_set
   else:
     exon_omit_list.append(1)
     #print gene_list

# Step 3: Go through the steps again to get the exonic sections. Each step
# becomes an 'exonic part'. The exonic part is associated with an
# aggregate gene, i.e., a gene set as determined in the previous step, 
# and a transcript set, containing all transcripts that occur in the step.
# The results are stored in the dict 'aggregates', which contains, for each
# aggregate ID, a list of all its exonic_part features.

aggregates = collections.defaultdict( lambda: list() )
i = -1
for iv, s in exons.steps( ):
   # Skip omitted steps
   i += 1
   if exon_omit_list[i]:
      continue
   # Take one of the gene IDs, find the others via gene sets, and
   # form the aggregate ID from all of them   
   gene_id = list(s)[0][0]
   assert set( gene_id for gene_id, transcript_id in s ) <= gene_sets[ gene_id ] 
   aggregate_id = '+'.join( gene_sets[ gene_id ] )
   # Make the feature and store it in 'aggregates'
   f = HTSeq.GenomicFeature( aggregate_id, "exonic_part", iv )   
   f.source = os.path.basename( sys.argv[1] )
   f.attr = {}
   f.attr[ 'gene_id' ] = aggregate_id
   transcript_set = set( ( transcript_id for gene_id, transcript_id in s ) )
   f.attr[ 'transcripts' ] = '+'.join( transcript_set )
   aggregates[ aggregate_id ].append( f )

# Step 4: For each aggregate, number the exonic parts

aggregate_features = []
for l in aggregates.values():
   name = l[0].name
   chrom = l[0].iv.chrom
   strand = l[0].iv.strand
   start = l[0].iv.start
   end = l[0].iv.end
   for i in xrange( 1, len(l)-1 ):
      assert l[i].name == name, str(l[i]) + str(l[0]) + " has wrong name"
      if l[i].iv.chrom != chrom:
        continue
        #raise ValueError, "Same name found on two chromosomes: %s, %s" % ( str(l[i]), str(l[i+1]) )
      if l[i].iv.strand != strand:
        continue
        #raise ValueError, "Same name found on two strands: %s, %s" % ( str(l[i]), str(l[i+1]) )
      assert l[i-1].iv.end <= l[i].iv.start, str(l[i-1]) + str(l[i]) + " starts too early"
      if start > l[i].iv.start :
        start = l[i].iv.start
      if end < l[i].iv.end :
        end = l[i].iv.end
   aggr_feat = HTSeq.GenomicFeature( name, "aggregate_gene", HTSeq.GenomicInterval( chrom, start, end, strand ) )
   aggr_feat.source = os.path.basename( sys.argv[1] )
   aggr_feat.attr = { 'gene_id': aggr_feat.name }
   for i in xrange( len(l) ):
      l[i].attr['exonic_part_number'] = "%03d" % ( i+1 )
   aggregate_features.append( aggr_feat )
      
      
# Step 5: Sort the aggregates, then write everything out

aggregate_features.sort( key = lambda f: ( f.iv.chrom, f.iv.start ) )

fout = open( out_file, "w" ) 
for aggr_feat in aggregate_features:
   fout.write( aggr_feat.get_gff_line() )
   for f in aggregates[ aggr_feat.name ]:
      if f.iv.strand != aggr_feat.iv.strand :
        continue
      if f.iv.chrom != aggr_feat.iv.chrom :
        continue 
      fout.write( f.get_gff_line() )

fout.close() 