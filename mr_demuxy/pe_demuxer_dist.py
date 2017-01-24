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



from itertools import izip
from biopython import FastqGeneralIterator
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from time import sleep

start_time = time.time()


## This does the command line parsing ##
## Probably should include an example input ## 

parser = argparse.ArgumentParser(add_help=False, description='This program splits \
combinatorially multiplexed paired-end seqs and outputs a folder for the R1 reads and \
another for the R2 reads. The two folders contain one fastq file per sample id.',
 usage='\npe_demuxer\
 -r1 r1_seqs.fastq -r2 r2_seqs.fastq -r1_bc r1_bc_file.txt -r2_bc -r2_bc_file.txt') #,add_help=False 
req = parser.add_argument_group('required arguments')
opt = parser.add_argument_group('optional arguments')

req.add_argument('-r1_bc',  
    help = 'BC file for R1 reads', 
    metavar = '',
    required=True)
req.add_argument('-r2_bc', 
    help = 'BC file for R2 reads',
    metavar = '' ,
    required=True)
req.add_argument('-r1', 
    help = 'R1 fastq file',
    metavar = '',
    required=True)
req.add_argument('-r2', 
    help = 'R2 fastq file',
    metavar = '',
    required=True)
#parser.add_argument('-h', '--help', action='store_true')
opt.add_argument("-h", "--help", action="help", help="show this help message and exit")
opt.add_argument('-o',
    help = 'name of output file or dir',
    default = 'pe_demuxer_output',
    metavar = '')
opt.add_argument('--batch_length',
    type=int, 
    default = 1500000,
    help='amount of seqs to be put in a batch \
    only used when check_headers is passed (Default: 1500000)',
    metavar='')
opt.add_argument('--forward_primer_len',
    help = 'base number of primer/adapter on R1 seq (default: 0)',
    type=int, 
    default=0,
    metavar = '')

opt.add_argument('--reverse_primer_len',
    help = 'base number of primer/adapter on R2 seq (default: 0)',
    type = int, 
    default = 0,
    metavar = '')
opt.add_argument('--keep_original_headers', 
    help = 'pass this arg to keep original headers in demultiplexed fasta/fastq file', 
    action='store_true')
opt.add_argument('--min_seq_length',
    help = 'minimum length of seq to be retained (default: 20)',
    metavar = '',
    default = 20)
opt.add_argument('--check_headers', 
    help = 'ensures the R1 and R2 read match before demultiplexing WARNING: this is super slow',
    action = 'store_true')


args = parser.parse_args()

example_usage = ('\n\t\t\tHow to use:\npe_demuxer\
 -r1 r1_seqs.fastq -r2 r2_seqs.fastq -r1_bc r1_bc_file.txt -r2_bc -r2_bc_file.txt')
example_r1_bc = ('\n\nThe -r1_bc argument specifies the barcode file corresponding to the r1 fastq file.\
The first column needs to correspond to letter (or whatever the front bc specifies)\
with the second being the barcode sequence.\
\nHere is an example of how this file should be formatted:' + \
'\n\tA GGTAC \n' + '\tB CAACAC \n' + '\tC ATCGGTT \n' + '---Note: the easiest way to make this file is\
to use excel and save it as a tab-delineated text file---\
\nAn example is provided by the example_data/forward_bcs.txt file.\n')
example_r2_bc = ('\nThe -r2_bc argument specifies the barcode file corresponding to the r2 fastq file.\
it should be formatted like the r1_bc file.\n')
example_r1_file = ('\nThe -r1 argument specifies the R1 seqs file.\
 See the example provided in the example_data dir for more information\n')
example_r2_file = ('\nThe -r2 argument specifies the R1 seqs file.\
 See the example provided in the example_data dir for more information\n')
example_keep = ('\n\nI dont actually think the --keep_original_headers\
does anything in this script, but it was too much trouble to remove.\n')
example_min = ('\nThe -min argument specifies what length the seq needs to be so that it is\
retained.\nDefault: 20')
example_check_header = ('\n\nThe --check_headers can handle situations where the seqs in the R1\
and R2 file are out of order and with reads missing from file.\nWARNING: this is going to be slow\
theres not much I can do about that, so I would recommend having your seqs in the same order.')
example_fp = '\n\nThe -fp argument specifies the how many bases after the barcode should be removed.\
The point of this is to remove primer/adapter sequences from reads without having to use any\
pattern matching. Thus, if you know how long the front primer/adapter sequences, you can just remove that\
many bases from the front of the read. If no value is passed, it defaults to 0.\n'
example_bp = '\n\nThe -bp argument is the same as the -fp, except it truncates so many bases from the end.\
Again, the point of this is to remove primer/adapter sequences without doing any sort of matching. This one \
also defaults to 0.\n'
example_output = ('\n\nThe -o just specifies the top dirs name. Two dirs within this dir\
 will be created, one for the R1 reads, and the other for the R2 reads.')



########################################################
''' These are the commands to make it go '''
########################################################

class PEDemux:
    def __init__(self):
        '''
        The PEDemux class takes care of Paired-end seqs.
        It demultiplexes by getting the front bits of off 
        the R1 and R2 read. If the reads are in the same 
        order and there is an R1 read for every R2 read, 
        the alt_pe_loop can be used. Else, use the slow_pe_loop
        which will only get the overlapping set of reads determined
        by header. This is unfortunately much slower.

        '''
        self.args = args
        self.util_functions_dist = util_functions_dist
        self.r1_bc_dict, self.r1_len = self.util_functions_dist.bc_dict_maker(self.args.r1_bc, False)
        self.r2_bc_dict, self.r2_len = self.util_functions_dist.bc_dict_maker(self.args.r2_bc, False)
        self.sam_ind = self.util_functions_dist.sample_index_maker(self.r1_bc_dict, self.r2_bc_dict)
        self.seq_number = 0
        
        self.r1_file_dict = self.util_functions_dist.outfile_opener(
            self.args.o,
            self.sam_ind,
            dir_opt='R1')
        
        self.r2_file_dict = self.util_functions_dist.outfile_opener(
            self.args.o,
            self.sam_ind,
            dir_opt='R2')
    
    
    # This creates a dict w/ with the seqs and a set w/ #
    # header info #
    def read_set_maker(self, seq_batch):
        '''
        Returns a dict with the header as key and the seq_rec 
        tuple as item. It also returns a set of headers which
        for determining which overlap.
        '''
        read_dict = {}
        cur_set = set()
        for seq_rec in seq_batch:
            hd, seq, qual = seq_rec
            try:
                head_id, pe_crap = hd.split()
                cur_set.add(head_id)
                read_dict[head_id] = seq_rec
            except ValueError:
                cur_set.add(hd)
                read_dict[hd] = seq_rec
        return(read_dict, cur_set)
    
    # function which returns demultiplexed seqs #
    def pe_demultiplexed_seqrec_maker(self, r1_seq_rec, r2_seq_rec):
        '''
        This basically returns a tuple with the sample_id
        and both seq_recs. Fixed the nasty little code duplication
        from earlier. 

        Possible stuff to be added in future:
        -fuzzy matching

        '''
        ## Borrowed Functions ##
        de_bc_loop = self.util_functions_dist.de_bc_loop
        demult_header = self.util_functions_dist.demult_header
                ##

        r1_head, r1_seq, r1_qual = r1_seq_rec 
        r2_head, r2_seq, r2_qual = r2_seq_rec
        r1_front_bit = r1_seq[0:self.r1_len]
        r2_front_bit = r2_seq[0:self.r2_len]
        
        r1_comb, r1_start = de_bc_loop(r1_front_bit, self.r1_bc_dict)
        r2_comb, r2_start = de_bc_loop(r2_front_bit, self.r2_bc_dict)
        
        cur_r1_seq = r1_seq[(r1_start + self.args.forward_primer_len):len(r1_seq)]
        cur_r1_qual = r1_qual[(r1_start + self.args.forward_primer_len):len(r1_seq)]
        cur_r2_seq = r2_seq[(r2_start + self.args.reverse_primer_len):len(r1_seq)]
        cur_r2_qual = r2_qual[(r2_start + self.args.reverse_primer_len):len(r1_seq)]

        sample_id, header = demult_header(r1_comb, r2_comb, self.seq_number, \
            r1_head, self.args.keep_original_headers)
           

        r1_demux_seqrec = r1_head, cur_r1_seq, cur_r1_qual
        r2_demux_seqrec = r2_head, cur_r2_seq, cur_r2_qual

        return (sample_id, r1_demux_seqrec, r2_demux_seqrec) 
   
   

    # Alternate demuxing loop #
    def alt_pe_loop(self, r1_batch, r2_batch):
        '''
        Basically does the normal demuxing with the assumption
        that the seqs in both files have a 1 to 1 correspondance
        and are in the same order. No intermediary steps are required
        as all 192 files are kept open.        

        '''
        ### Borrowed Functions ####
        de_bc_loop = self.util_functions_dist.de_bc_loop
        demult_header = self.util_functions_dist.demult_header
        sample_index_maker = self.util_functions_dist.sample_index_maker
                    ###
        
        for (r1_cur_rec, r2_cur_rec) in izip(r1_batch, r2_batch):
            self.seq_number = self.seq_number + 1
            if self.seq_number % 1000 == 0:
                sys.stdout.write('Demultiplexed Seqs: ' '{0}\r'.format(self.seq_number),) 
                sys.stdout.flush()     
            
            sample_id, r1_demux_seqrec, r2_demux_seqrec = self.pe_demultiplexed_seqrec_maker(
                r1_cur_rec, 
                r2_cur_rec)
            
            r1_head, cur_r1_seq, cur_r1_qual = r1_demux_seqrec
            r2_head, cur_r2_seq, cur_r2_qual = r2_demux_seqrec

            if (len(cur_r1_seq) < self.args.min_seq_length) or (len(cur_r2_seq) < self.args.min_seq_length):
                continue
            else:
                
                r1_cur_fastq = self.r1_file_dict[sample_id]
                r1_cur_fastq.write('@' + r1_head + \
                    '\n' + cur_r1_seq + '\n' + '+' +'\n' + cur_r1_qual + '\n')
                r2_cur_fastq = self.r2_file_dict[sample_id]
                r2_cur_fastq.write('@' + r2_head + \
                    '\n' + cur_r2_seq + '\n' + '+' +'\n' + cur_r2_qual + '\n')


    # main demuxing loop
    def slow_pe_loop(self, r1_seqdict, r1_set, r2_seqdict, r2_set):
        '''
        Basically does the same as the above function, but uses a set of 
        headers from both files so that only reads with headers in both
        will be used.
        
        Really need to fix the duplication problem

        '''
        ### Borrowed Functions ### 
        sample_index_maker = self.util_functions_dist.sample_index_maker
                ###
        
        pe_set = r1_set & r2_set
        pe_list = list(pe_set)
        r1_reads_dict = sample_index_maker(self.r1_bc_dict, self.r2_bc_dict)
        r2_reads_dict = sample_index_maker(self.r1_bc_dict, self.r2_bc_dict)
        
        for seq_id in pe_list:
            self.seq_number = self.seq_number + 1
            if self.seq_number % 1000 == 0:
                sys.stdout.write('Demultiplexed Seqs: ' '{0}\r'.format(self.seq_number),) 
                sys.stdout.flush()

            r1_cur_rec = r1_seqdict[seq_id]
            r2_cur_rec = r2_seqdict[seq_id]
            
            sample_id, r1_demux_seqrec, r2_demux_seqrec =\
            pe_demultiplexed_seqrec_maker(r1_cur_rec, r2_cur_rec)
            
            r1_reads_dict[sample_id].append(r1_proc_seq)
            r2_reads_dict[sample_id].append(r2_proc_seq)

        return(r1_reads_dict, r2_reads_dict)


    def main_loop(self):
        ## Borrowed functions ##
        ind_fastq_writer = self.util_functions_dist.ind_fastq_writer
                        ##

        r1_iter = FastqGeneralIterator(open(self.args.r1, 'rU'))
        r2_iter = FastqGeneralIterator(open(self.args.r2, 'rU'))

        if self.args.check_headers is True:
            r1_generator = self.util_functions_dist.batch_iterator(
                r1_iter, 
                self.args.batch_length)
            r2_generator = self.util_functions_dist.batch_iterator(
                r2_iter, 
                self.args.batch_length)
            for (r1_batch, r2_batch) in izip(r1_generator, r2_generator):
                cur_r1_dict, cur_r1_set = self.read_set_maker(r1_batch)
                cur_r2_dict, cur_r2_set = self.read_set_maker(r2_batch)
                r1_demux_dict, r2_demux_dict = self.slow_pe_loop(
                    cur_r1_dict, 
                    cur_r1_set, 
                    cur_r2_dict, 
                    cur_r2_set)
                
                ind_fastq_writer(r1_demux_dict, self.r1_path, 'R1')
                ind_fastq_writer(r2_demux_dict, self.r2_path, 'R2')



        if self.args.check_headers is False:
            self.alt_pe_loop(
                r1_iter,
                r2_iter
                )
            for key, r1_file in self.r1_file_dict.iteritems():
                r1_file.close()
            
            for key, r2_file in self.r2_file_dict.iteritems():
                r2_file.close()


        

#if __name__ == '__main__':
#    pe_class = PEDemux()
#    pe_class.main_loop()
print('\n' "Program Runtime: " + str(round(((time.time() - start_time))/60, 3)) + ' (min)' + '\n')


    





            
                
    
        
