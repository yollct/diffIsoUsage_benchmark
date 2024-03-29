#!/usr/bin/env python2
# encoding: utf-8
"""
cuffmerge.py
Created by Cole Trapnell on 2011-03-17.
Copyright (c) 2011 Cole Trapnell. All rights reserved.
"""

import sys
import getopt
from datetime import datetime, date, time
import shutil
import subprocess
import errno
import os
import tempfile
import warnings
import types

help_message = '''
cuffmerge takes two or more Cufflinks GTF files and merges them into a
single unified transcript catalog.  Optionally, you can provide the script
with a reference GTF, and the script will use it to attach gene names and other
metadata to the merged catalog.
Usage:
    cuffmerge [Options] <assembly_GTF_list.txt>
Options:
    -h/--help                               Prints the help message and exits
    -o                     <output_dir>     Directory where merged assembly will be written  [ default: ./merged_asm  ]
    -g/--ref-gtf                            An optional "reference" annotation GTF.
    -s/--ref-sequence      <seq_dir>/<seq_fasta> Genomic DNA sequences for the reference.
    --min-isoform-fraction <0-1.0>          Discard isoforms with abundance below this       [ default:           0.05 ]
    -p/--num-threads       <int>            Use this many threads to merge assemblies.       [ default:             1  ]
    --keep-tmp                              Keep all intermediate files during merge
'''

output_dir = "./merged_asm/"
logging_dir = output_dir + "logs/"
run_log = None
run_cmd = None
tmp_dir = output_dir + "/tmp/"
bin_dir = sys.path[0] + "/"
run_meta_assembly = True
fail_str = "\t[FAILED]\n"
params = None

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class TestParams:

    class SystemParams:
        def __init__(self,
                     threads,
                     keep_tmp):
            self.threads = threads
            self.keep_tmp = keep_tmp

        def parse_options(self, opts):
            for option, value in opts:
                if option in ("-p", "--num-threads"):
                    self.threads = int(value)
                if option in ("--keep-tmp"):
                    self.keep_tmp = True

        def check(self):
            pass

    def __init__(self):
        self.system_params = self.SystemParams(1,               # threads
                                               False)           # keep_tmp
        self.ref_gtf = None
        self.fasta = None
        self.min_isoform_frac = 0.05

    def check(self):
        self.system_params.check()

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:],
                                       "hvp:o:g:M:s:q:F:",
                                       ["version",
                                        "help",
                                        "ref-sequence=",
                                        "ref-gtf=",
                                        "output-dir=",
                                        "num-threads=",
                                        "keep-tmp",
                                        "min-isoform-fraction="])
        except getopt.error, msg:
            raise Usage(msg)

        self.system_params.parse_options(opts)

        global output_dir
        global logging_dir
        global tmp_dir

        # option processing
        for option, value in opts:
            if option in ("-v", "--version"):
                print "merge_cuff_asms v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-g", "--ref-gtf"):
                self.ref_gtf = value
            if option in ("-s", "--ref-sequence"):
                self.fasta = value
            if option in ("-F", "--min-isoform-fraction"):
                self.min_isoform_frac = float(value)
            if option in ("-o", "--output-dir"):
                output_dir = value + "/"
                logging_dir = output_dir + "logs/"
                tmp_dir = output_dir + "tmp/"

        return args


def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def prepare_output_dir():

    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:
        os.makedirs(output_dir)

    #print >> sys.stderr, "Checking for %s", logging_dir
    if os.path.exists(logging_dir):
        pass
    else:
        #print >> sys.stderr, "Creating %s", logging_dir
        os.makedirs(logging_dir)

    if os.path.exists(tmp_dir):
        pass
    else:
        os.makedirs(tmp_dir)

def formatTD(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)

def tmp_name(prefix):
    tmp_root = output_dir + "tmp/"
    if os.path.exists(tmp_root):
        pass
    else:
        os.mkdir(tmp_root)
    return tmp_root + prefix + os.tmpnam().split('/')[-1]

def cufflinks(out_dir,
              sam_file,
              min_isoform_frac,
              gtf_file=None,
              extra_opts=["-q", "--overhang-tolerance", "200", "--library-type=transfrags",  "-A","0.0", "--min-frags-per-transfrag", "0", "--no-5-extend"],
              lsf=False,
              curr_queue=None):
    if gtf_file != None:
        print >> sys.stderr, "[%s] Quantitating transcripts" % (right_now())
    else:
        print >> sys.stderr, "[%s] Assembling transcripts" % (right_now())

    cmd = ["cufflinks"]

    if out_dir != None and out_dir != "":
        cmd.extend(["-o", out_dir])

    cmd.extend(["-F", str(min_isoform_frac)])

    if gtf_file != None:
        cmd.extend(["-g", gtf_file])

    if extra_opts != None:
        cmd.extend(extra_opts)
    global params
    # Run Cufflinks with more than one thread?
    cmd.extend(["-p", str(params.system_params.threads)])

    cmd.append(sam_file)

    try:
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cufflinks"
            exit(1)
    # cufflinks not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cufflinks not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def cuffcompare(prefix, ref_gtf, fasta, cuff_gtf):

    print >> sys.stderr, "[%s] Comparing reference %s to assembly %s" % (right_now(), ref_gtf, cuff_gtf)
    cmd = ["cuffcompare"]

    if  prefix != None:
        cmd.extend(["-o", prefix])
    if  ref_gtf != None:
        cmd.extend(["-r", ref_gtf])
    if  fasta != None:
        cmd.extend(["-s", fasta])
    if type(cuff_gtf) == types.ListType:
        for g in cuff_gtf:
            cmd.extend([g])
    else:
        cmd.extend([cuff_gtf])

    try:
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def gtf_to_sam(gtf_filename):

    sam_out = tmp_name("gtf2sam_")

    cmd = ["/nfs/home/students/chit/gtf_to_sam"]
    cmd.append("-F")
    cmd.append(gtf_filename)
    cmd.append(sam_out)
    try:
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute gtf_to_sam"
            exit(1)
    # gtf_to_sam not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: gtf_to_sam not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
    return sam_out

def test_input_files(filename_list):
    """This function takes a file that contains a list of GTF files,
       tests accessibility of each, and returns the list of filenames"""

    OK = True
    input_files = []
    for line in filename_list:
        line = line.strip()

        # Skip comment line
        if len(line) == 0 or line[0] == "#":
            continue
        try:
            g = open(line,"r")
            input_files.append(line)

        except OSError, o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                print >> sys.stderr, fail_str, "Error: could not open %s" % line
            OK = False
    if not OK:
        sys.exit(1)
    return input_files

def convert_gtf_to_sam(gtf_filename_list):
    """This function takes a list of GTF files, converts them all to
       temporary SAM files, and returns the list of temporary file names."""
    print >> sys.stderr, "[%s] Converting GTF files to SAM" % (right_now())
    OK = True
    sam_input_filenames = []
    for line in gtf_filename_list:
        try:
            sam_out = gtf_to_sam(line)
            sam_input_filenames.append(sam_out)
        except OSError, o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                print >> sys.stderr, fail_str, "Error: could not open %s" % line
            OK = False
    if not OK:
        sys.exit(1)
    return sam_input_filenames

def merge_sam_inputs(sam_input_list, header):
    sorted_map_name = tmp_name( "mergeSam_")

    sorted_map = open(sorted_map_name, "w")

    #print header

    # The header was built from a dict keyed by chrom, so
    # the records will be lexicographically ordered and
    # should match the BAM after the sort below.
    print >> sorted_map, header,

    sorted_map.close()
    sort_cmd =["sort",
               "-k",
               "3,3",
               "-k",
               "4,4n",
               "--temporary-directory="+tmp_dir]
    sort_cmd.extend(sam_input_list)

    print >> run_log, " ".join(sort_cmd), ">", sorted_map_name
    subprocess.call(sort_cmd,
                   stdout=open(sorted_map_name, "a"))
    return sorted_map_name

def compare_to_reference(meta_asm_gtf, ref_gtf, fasta):
    global tmp_dir
    print >> sys.stderr, "[%s] Comparing against reference file %s" % (right_now(), ref_gtf)
    ref_str = ""
    if ref_gtf != None:
        ref_str = " -r %s " % ref_gtf

    if fasta != None:
        comp_cmd = '''cuffcompare -o %s -C -G %s -s %s %s''' % (tmp_dir, ref_str, fasta, meta_asm_gtf)
    else:
        comp_cmd = '''cuffcompare -o %s -C -G %s %s''' % (tmp_dir, ref_str, meta_asm_gtf)

    #cmd = bsub_cmd(comp_cmd, "/gencode_cmp", True, job_mem=8)
    cmd = comp_cmd

    try:
        print >> run_log, cmd
        ret = subprocess.call(cmd,shell=True)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
        #tmap_out = meta_asm_gtf.split("/")[-1] + ".tmap"
        tfpath, tfname = os.path.split(meta_asm_gtf)
        if tfpath: tfpath+='/'
        tmap_out = tfpath+'.'+tfname+".tmap"
        return tmap_out
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def select_gtf(gtf_in_filename, ids, gtf_out_filename):
    f_gtf = open(gtf_in_filename)
    #print >> sys.stderr, "Select GTF: Ids are: "
    #print >> sys.stderr, ids
    #print >> sys.stderr, "reference gtf file name:"
    #print >> sys.stderr, gtf_in_filename
    out_gtf = open(gtf_out_filename, "w")
    for line in f_gtf:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        attrs = cols[8]
        attr_cols = attrs.split(';')
        for col in attr_cols:
            if col.find("transcript_id") != -1:
                first_quote = col.find('"')
                last_quote = col.find('"', first_quote + 1)
                transcript = col[first_quote + 1:last_quote]
                #print >> sys.stderr, transcript
                if transcript in ids:
                    print >> out_gtf, line


def merge_gtfs(gtf_filenames, merged_gtf, ref_gtf=None):
    print >> sys.stderr, "[%s] Merging linc gtf files with cuffcompare" % (right_now())
    cmd = ["cuffcompare"]

    cmd.extend(["-o", merged_gtf])
    if ref_gtf != None:
        cmd.extend(["-r", ref_gtf])

    cmd.extend(gtf_filenames)
    cmd = " ".join(cmd)
    #cmd = bsub_cmd(cmd, "/merge_gtf", True, job_mem=8)

    try:
        print >> run_log, cmd
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
        return merged_gtf + ".combined.gtf"
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def compare_meta_asm_against_ref(ref_gtf, fasta_file, gtf_input_file, class_codes=["c", "i", "r", "p", "e"]):
    #print >> sys.stderr, "Cuffcmpare all assemblies GTFs"

    tmap = compare_to_reference(gtf_input_file, ref_gtf, fasta_file)

    #print >> sys.stderr, "Cuffcmpare all assemblies GTFs : filter %s" % ",".join(class_codes)
    selected_ids= set([])
    f_tmap = open(tmap)
    #out = open("tmp_meta_asm_selectedIds.txt", "w")
    for line in f_tmap:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 5:
            continue
        class_code = cols[2]
        name = cols[4]
        if class_code not in class_codes:
            selected_ids.add(name)

    global output_dir
    asm_dir = output_dir

    if os.path.exists(asm_dir):
        pass
    else:
        os.mkdir(asm_dir)
    current_asm_gtf = output_dir +"transcripts.gtf"
    select_gtf(current_asm_gtf, selected_ids, output_dir + "/merged.gtf")
    mtmap = compare_to_reference(output_dir + "/merged.gtf", ref_gtf, fasta_file)
    if os.path.exists(mtmap):
        os.remove(mtmap)
    if os.path.exists(mtmap.split(".tmap")[0]+".refmap"):
        os.remove(mtmap.split(".tmap")[0]+".refmap")

    global tmp_dir
    shutil.move(tmp_dir+".combined.gtf", output_dir + "/merged.gtf")

#    os.remove("tmp_meta_asm.combined.gtf")
    if os.path.exists(tmap):
        os.remove(tmap)
    if os.path.exists(tmap.split(".tmap")[0]+".refmap"):
        os.remove(tmap.split(".tmap")[0]+".refmap")
    #tmp_dir = asm_dir
    #tmp_files = os.listdir(tmp_dir)
    #for t in tmp_files:
    #    os.remove(tmp_dir+t)
    #os.rmdir(tmp_dir)

#os.remove("tmp_meta_asm.tmap")

def get_version():
    return "1.0.0"

def get_gtf_chrom_info(gtf_filename, known_chrom_info=None):
    gtf_file = open(gtf_filename)
    if known_chrom_info == None:
        chroms = {}
    else:
        chroms = known_chrom_info
    for e, line in enumerate(gtf_file):
        line = line.strip()

        if len(line) == 0 or line[0] == "#":
            continue
        cols = line.split('\t')
        if len(cols) < 8:
            continue
        chrom = cols[0]
        left = int(cols[3])
        if cols[4]=="Cufflinks":
            print(gtf_filename, e)
            continue 
        right = int(cols[4])
        bounds = chroms.setdefault(chrom, [9999999999,-1])
        if bounds[0] > left:
            bounds[0] = left
        if bounds[1] < right:
            bounds[1] = right
    return chroms

def header_for_chrom_info(chrom_info):
    header_strs = ["""@HD\tVN:1.0\tSO:coordinate"""]
    chrom_list = [(chrom, limits) for chrom, limits in chrom_info.iteritems()]
    chrom_list.sort(lambda x,y: cmp(x[0],y[0]))
    #print chrom_list
    for chrom, limits in chrom_list:
        line = "@SQ\tSN:%s\tLN:\t%d" % (chrom, limits[1])
        header_strs.append(line)
    header_strs.append("@PG\tID:cuffmerge\tVN:1.0.0\n")
    header = "\n".join(header_strs)
    return header

def main(argv=None):
    
    # Set this so sort is consistent across platforms and python versions 
    # (and always character-lexicographic)
    os.environ['LC_ALL']='C'
    warnings.filterwarnings("ignore", "tmpnam is a potential security risk")
    global params
    params = TestParams()

    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()

        if len(args) < 1:
            raise(Usage(help_message))

        global run_log
        global run_cmd

        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning transcriptome assembly merge" % (right_now())
        print >> sys.stderr, "-------------------------------------------"
        print >> sys.stderr

        start_time = datetime.now()
        prepare_output_dir()

        run_log = open(logging_dir + "run.log", "w", 0)
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd

        transfrag_list_file = open(args[0], "r")

        if params.ref_gtf != None:
            test_input_files([params.ref_gtf])
        else:
            print >> sys.stderr, "Warning: no reference GTF provided!"

        # Check that all the primary assemblies are accessible before starting the time consuming stuff
        gtf_input_files = test_input_files(transfrag_list_file)

        all_gtfs = []
        all_gtfs.extend(gtf_input_files)
        if params.ref_gtf != None:
            all_gtfs.append(params.ref_gtf)
        chrom_info = {}
        for gtf in all_gtfs:
            chrom_info = get_gtf_chrom_info(gtf, chrom_info)


        header = header_for_chrom_info(chrom_info)

        #Meta assembly option:
        global run_meta_assembly
        if run_meta_assembly:
            # Convert the primary assemblies to SAM format
            sam_input_files = convert_gtf_to_sam(gtf_input_files)
            # Merge the primary assembly SAMs into a single input SAM file
            merged_sam_filename = merge_sam_inputs(sam_input_files, header)
            # Run cufflinks on the primary assembly transfrags to generate a meta-assembly
            cufflinks(output_dir, merged_sam_filename, params.min_isoform_frac, params.ref_gtf)
            compare_meta_asm_against_ref(params.ref_gtf, params.fasta, output_dir+"/transcripts.gtf")
        #Meta Cuffcompare option:
        else:
            cuffcompare_all_assemblies(gtf_input_files) #FIXME: where is this function ?


        if not params.system_params.keep_tmp:
            tmp_files = os.listdir(tmp_dir)
            for t in tmp_files:
                os.remove(tmp_dir+t)
            os.rmdir(tmp_dir)
            os.remove(output_dir + "/transcripts.gtf")
            os.remove(output_dir + "/skipped.gtf")
            os.remove(output_dir + "/genes.fpkm_tracking")
            os.remove(output_dir + "/isoforms.fpkm_tracking")
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2


if __name__ == "__main__":
    sys.exit(main())