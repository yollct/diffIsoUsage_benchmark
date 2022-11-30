library("optparse")

option_list = list(
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-g", "--ngenes"), type="numeric", default="1000")
); 

suppressMessages(library(ASimulatoR))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

gtf_file = system.file('extdata', 'Homo_sapiens.GRCh38.99.21.gtf', package = 'ASimulatoR')
exon_superset = get_exon_supersets(gtf_file)

# define your input_dir, where the annotation gtf (or the exon supersets if you have already created them) and the genome fasta files (one per chromosome) are located
# here we will use the example data
input_dir = system.file('extdata', package = 'ASimulatoR')

# define, how many groups and samples per group you analyze. Here we create a small experiment with two groups with one sample per group:
num_reps = c(3,3)

# define your outdir with NO slash
outdir = opt$out

# in this example we use relative frequencies
# here we produce eight variants with one of each AS events as well as one variant containing every event
# if probs_as_freq was FALSE, a random number would be drawn for each event-superset combination and only if it was smaller than 1/9 the AS event would be created
probs_as_freq = F
event_freq = 
  setNames(rep(1/9, 9),
           c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee', 'es,ir,mes,a3,a5,afe,ale,mee'))
as_events = c('es','ir','a3','a5','mes','mee','ale','afe')
event_probs = rep(1/(length(as_events) + 1), length(as_events))
names(event_probs) = as_events

max_genes = 500 #number of genes. NULL - use all compatible exon supersets
seq_depth = 2e06 #sequencing depth i.e. the number of reads per sample
multi_events_per_exon = T
error_rate = 0.001
readlen = 76
seed=142
ncores=1

# we use the previously created superset to simulate splice variants from, since it is saved in the same directory as the gtf
# if no superset is found, a new one will be created
params = list(
  seed = seed,
  ncores = ncores,
  input_dir = input_dir,
  event_probs = event_probs,
  outdir = outdir,
  seq_depth = seq_depth,
  max_genes = max_genes,
  error_rate = error_rate,
  readlen = readlen,
  multi_events_per_exon = multi_events_per_exon,
  probs_as_freq = probs_as_freq,
  num_reps = num_reps
)


### run simulator ----
library(ASimulatoR)
do.call(simulate_alternative_splicing, params) #simulate RNA-Seq datasets








