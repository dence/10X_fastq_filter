#Daniel Ence
#June 27, 2018

import re
import sys
import os.path
import gzip
from gzip import open as gzopen
from Bio import SeqIO, bgzf
import multiprocessing
from itertools import izip
import argparse

parser = argparse.ArgumentParser(description="Filter fastq files for 10X R2 repeat seqs")

parser.add_argument("-1", "--read1", help="Paired-end fastq file one")
parser.add_argument("-2", "--read2", help="Paired-end fastq file two")
parser.add_argument("-p","--num-threads", help="Number of threads to start",type=int)
parser.add_argument("-o", "--output", help="Output filename (interpreted as basename for PE)")
parser.add_argument("-f", "--filter_file",help="File containing sequences to filter for ")
parser.add_argument("--gzip", help="gzip compress the output",action="store_true")

args = parser.parse_args()

filename, file_extension = os.path.splitext(args.read1)
if file_extension == ".gz":
    infile_read1 = gzip.open(args.read1,"rb")
elif file_extension == ".fastq":
    infile_read1 = open(args.read1,"r")
else:
    sys.exit("Fatal: unknown file extension for filename %" % args.read1)

file_extension, file_extension = os.path.splitext(args.read2)
if file_extension == ".gz":
    infile_read2 = gzip.open(args.read2,"rb")
elif file_extension == ".fastq":
    infile_read2 = open(args.read2,"r")
else:
    sys.exit("Fatal: uknown file extension for filename %" % args.read2)

if args.gzip:
    outfile_read1 = gzip.open(args.output + ".1.fastq.gz","wb")
    outfile_read2 = gzip.open(args.output + ".2.fastq.gz","wb")
else:
    outfile_read1 = open(args.output + ".1.fastq", "w")
    outfile_read2 = open(args.output + ".2.fastq", "w")

sys.stderr.write("Filtering reads from paired files %s and %s" % (args.read1,args.read2))

filter_list = [i.strip() for i in open(args.filter_file)]

sys.stderr.write("Read %d bad sequences to filter.\n" % len(filter_list))

def filter_pair(lines):
    #lines is a set of four lines from both files that are stitched together. Order is:
    # file1line1, file2line2, file1line2, file2line2, etc...

    seqid1 = lines[0].split()[0]
    seqid2 = lines[1].split()[0]

    sequence2 = lines[2].split()[0]

    print "sequence2 is %s" % sequence2

    if seqid1 != seqid2:
        print "Warning: out of sync pairs detected with ids %s %s" % (seqid1,seqid2)

    for badSeq in filter_list:
        if re.search(badSeq,sequence2):
            print "Found this:\t" + badSeq
            print "in this:\t" + sequence2
            return False
    return lines

n_reads = 0
n_filtered = 0

pool = multiprocessing.Pool(processes=args.num_threads)

res = pool.imap(filter_pair, izip(*[infile_read1, infile_read2]*4))

while 1:
    if not n_reads % 1000000:
        sys.stderr.write("Processed %d reads.\n" % n_reads)
    try:
        lines = res.next()

        n_reads += 1
        if lines:
            outfile_read1.write(lines[0] + lines[2] + lines[4] + lines[6])
            outfile_read2.write(lines[1] + lines[3] + lines[5] + lines[7])

        else:
            n_filtered += 1

    except StopIteration:
        break

pool.close()
pool.join()

sys.stderr.write("Processed %d reads, filtered %d reads. (%.2f%%).\n" % (n_reads, n_filtered, n_filtered/float(n_reads)))
