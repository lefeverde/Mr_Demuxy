# -*- coding: utf-8 -*-
'''
Created on Fri Feb 12 19:00:29 2016

@author: Daniel Lefever
@email: lefeverd@uga.edu


This program was written using Biopython:

Cock, P.J.A. et al. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878

'''
#!/usr/bin/python


## These things make the program go ##
from __future__ import division


import sys
import argparse
import time
import os
import gc
import tempfile
import util_functions_dist


from biopython import FastqGeneralIterator
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from time import sleep

start_time = time.time()


## This does the command line parsing ##
## Probably should include an example input ## 

parser = argparse.ArgumentParser(add_help=False, 
    description='This program demultiplexes combinatorially barcoded reads.',
    usage='\nmerged_demuxer -f forward_bcs.txt -r reverse_bcs.txt \
-i example_merged.fastq  -o demultiplexed_seqs.fasta') #, add_help=False
req = parser.add_argument_group('required arguments')
opt = parser.add_argument_group('optional arguments')
req.add_argument('-f',
     help='forward bcs file',
     metavar = '',
     required=True) 
req.add_argument('-r', 
    help = 'reverse bcs file',
    metavar = '',
    required=True) 
req.add_argument('-i', 
    help = 'input merged fastq file',
    metavar = '',
    required=True) 
opt.add_argument('-o',
    help = 'name of output file or folder (default: merged_demuxer_out)',
    default = 'merged_demuxer_out',
    metavar = '') 
opt.add_argument('--forward_primer_len',
    help = 'base number of forward primer/adapter (default: 0)',
    type=int, 
    default=0,
    metavar = '')
opt.add_argument('--reverse_primer_len',
    help = 'base number of reverse primer/adapter (default: 0)',
    type = int, 
    default = 0,
    metavar = '')

opt.add_argument('--no_rev_comp_reverse_bcs',
    action = 'store_false',
    help ='prevents rev bcs from reverse complementing'
    )
opt.add_argument("-h", "--help", action="help", help="show this help message and exit")

opt.add_argument('--keep_original_headers', 
    help = 'pass this arg to keep original headers in demultiplexed fasta/fastq file', 
    action='store_true')
opt.add_argument('--write_qual', 
    action='store_true',
    help='write a qual file containing quality scores')
opt.add_argument('--individual_fastq', 
    action='store_true',
    help='demultiplex by writing a fastq file for each sample ID')
opt.add_argument('--min_seq_length',
    help = 'minimum length of seq to be retained (default: 20)',
    metavar = '',
    default = 20)
opt.add_argument('--file_extension', 
    default='fasta',
    help='format of the output file (default: fasta)',
    metavar = '')



example_usage = ('\n\t\t\tHere is an example of how to use this program:\
\n-f forward_bcs.txt -r reverse_bcs.txt \
-i example_merged.fastq  -o demultiplexed_seqs.fasta')    
example_forward_bcs = '\nThe -f argument specifies the forward barcode file.\
The first column needs to correspond to letter (or whatever the front bc specifies)\
with the second being the barcode sequence.\
\nHere is an example of how this file should be formatted:' + \
'\n\tA GGTAC \n' + '\tB CAACAC \n' + '\tC ATCGGTT \n' + '---Note: the easiest way to make this file is\
to use excel and save it as a tab-delineated text file---\
\nAn example is provided by the example_data/forward_bcs.txt file.\n'
example_reverse_bcs = '\n\nThe -r argument specifies the reverse barcode file. like the forward_bcs file, \
it needs to have the first column to correspond to number (or whatever the back bc specifies) \
with the second corresponding to barcode sequence\
\nHere is example of how this one should be formatted:' + \
'\n1 AGGAA \n' + '2 GAGTGG\n' + '3 CCACGTC\n \
An example is provided by the example_data/reverse_bcs.txt file\n'
example_input_fasta = '\n\nThe -i argument specifies the input merged fastq file. This is the unmultiplexed file \
created by merging pair-ended fastq files. Note: this script only accepts merged fastq files.\
\nFor paired-end files, use the pe_demuxer script. \nThe merged fastq should looks somehting like this:\n\n\
@M02849:138:000000000-AJTC8:1....\nGGTACCCTACGG.....GGTAGTCTTCCT......'
example_output_fasta = '\n\nThe -o option just specifies what to the output file or dir will be called.\
 Call it whatever you want.\n'
example_fp = '\nThe -fp argument specifies the how many bases after the barcode should be removed.\
The point of this is to remove primer/adapter sequences from reads without having to use any\
pattern matching. Thus, if you know how long the front primer/adapter sequences, you can just remove that\
many bases from the front of the read. If no value is passed, it defaults to 0.\n'
example_bp = '\nThe -bp argument is the same as the -fp, except it truncates so many bases from the end.\
Again, the point of this is to remove primer/adapter sequences without doing any sort of matching. This one \
also defaults to 0.\n'
example_ext = '\nThe -ext argument specifies the output file type. The default setting is fasta. \
to get a fastq file, pass -ext fastq.'
example_ind = '\nThe -ind argument will split the fastq into individual fq files with one sample per file\n '
example_keep = '\nThe -keep argument tells the program whether to keep the orignal illumina fastq header \
information when a single fastq, fasta, or fasta and qual are the output.\
The new header becomes @<orignal_header_from_sequencers> <demultiplexed_seqID>'
example_write_qual = '\n\nThe --write_qual argument allows the quality scores of each seq to written into a \
seperate file. Simply pass the argument to use.'
args = parser.parse_args()







########################################################
''' These are the commands to make it go '''
########################################################

class MergedDemuxer:
    def __init__(self):
        '''
        Class for handling the demultiplexing of merged reads
        
        '''
        self.args = args
        self.util_functions = util_functions_dist
        self.fbc_dict, self.fbc_len = self.util_functions.bc_dict_maker(
            self.args.f,
            False)
        self.rbc_dict, self.rbc_len = self.util_functions.bc_dict_maker(
            self.args.r,
            self.args.no_rev_comp_reverse_bcs)
        self.sam_ind = self.util_functions.sample_index_maker(self.fbc_dict, self.rbc_dict)
        self.seq_number = 0
        self.infile_extension = self.util_functions.infile_extension(self.args.i)
        self.out_path = os.path.abspath(self.args.o)
       

    def demultiplexed_seqrec_maker(self, seq_rec):
        de_bc_loop = self.util_functions.de_bc_loop
        trimmed_read = self.util_functions.trimmed_read
        demult_header = self.util_functions.demult_header

        orig_header, seq, qual = seq_rec
        front_subseq = seq[0:self.fbc_len]
        back_subseq = seq[-self.rbc_len:]
        
        front_id, start_pos = de_bc_loop(front_subseq, self.fbc_dict)
        
        back_id, end_pos = de_bc_loop(back_subseq, self.rbc_dict)

        cur_seq, cur_qual = trimmed_read(seq, qual, start_pos, end_pos,\
            self.args.forward_primer_len, self.args.reverse_primer_len, self.infile_extension)
        
        sample_id, cur_header = demult_header(front_id, back_id, self.seq_number, orig_header,\
            self.args.keep_original_headers)

        demult_seq_rec = (sample_id, orig_header, cur_header, cur_seq, cur_qual)
        
        return demult_seq_rec


    def merged_demux_loop(self, batch, file_dict=None, out_file=None, out_qual=None):
        '''
        This is the main demuxing looping 

        '''
        demultiplexed_seqrec_maker = self.demultiplexed_seqrec_maker
       
        for seq_rec in batch:
            self.seq_number = self.seq_number + 1
            if self.seq_number % 1000 == 0:
                sys.stdout.write('Demuxed Seqs: ' '{0}\r'.format(self.seq_number),) 
                sys.stdout.flush()

            sample_id, orig_header, cur_header, cur_seq, cur_qual = demultiplexed_seqrec_maker(seq_rec)
            
            if self.args.individual_fastq is False:

                if self.args.file_extension == 'fasta':
                    if len(cur_seq) < self.args.min_seq_length:
                        continue
                    elif self.args.write_qual is False:
                        out_file.write('>' + cur_header + '\n' + cur_seq + '\n')
                        continue
                    else:
                        out_file.write('>' + cur_header + '\n' + cur_seq + '\n')
                        out_qual.write('>' + cur_header + '\n' + cur_qual + '\n')
                        continue
                       
                if self.args.file_extension == 'fastq':
                    if len(cur_seq) < self.args.min_seq_length:
                        continue
                    else:
                        out_file.write('@' + cur_header + '\n' + cur_seq \
                        + '\n' + '+' +'\n' + cur_qual + '\n')
                        continue

            if self.args.individual_fastq is True:
                cur_fastq = file_dict[sample_id]
                if len(cur_seq) < self.args.min_seq_length:
                    continue
                else:
                    cur_fastq.write('@' + orig_header + \
                    '\n' + cur_seq + '\n' + '+' +'\n' + cur_qual + '\n')
            
         
        return file_dict

    def main_loop(self):
        '''
        This is going demultiplex and split everything all in one go,
        since the files are broken into bits anyhow 
        '''
        fq_generator = FastqGeneralIterator(open(self.args.i))
        merged_demux_loop = self.merged_demux_loop

        
        if self.args.individual_fastq is True:
            file_dict = self.util_functions.outfile_opener(self.out_path, self.sam_ind)
            merged_demux_loop(
                fq_generator,
                file_dict = file_dict
                )
            for key, file in file_dict.iteritems():
                file.close()


        if self.args.individual_fastq is False:
            if self.args.write_qual is True:
                out_file_open = open((self.out_path + '.' + self.args.file_extension), 'w+')
                out_qual_open = open((self.out_path + '.qual'), 'w+')
                merged_demux_loop(
                    fq_generator, 
                    out_file=out_file_open, 
                    out_qual=out_qual_open
                    )
                out_qual_open.close()
                out_file_open.close()

        if self.args.individual_fastq is False:
            if self.args.write_qual is False:
                out_file_open = open((self.out_path + '.' + self.args.file_extension), 'w+')
                merged_demux_loop(
                    fq_generator, 
                    out_file=out_file_open
                    )         
                out_file_open.close() 
               

#if __name__ == '__main__':

#    sys.stdout.write('\n')
#    sys.stdout.flush()
#    merged_class = MergedDemuxer()
#    merged_class.main_loop()
#    sys.dont_write_bytecode = True

print('\n' "Program Runtime: " + str(round(((time.time() - start_time))/60, 3)) + ' (min)' + '\n')


    





            
                
    
        
        