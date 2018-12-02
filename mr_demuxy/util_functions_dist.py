#!/usr/bin/python

## These things make the program go ##


import sys
import os
import itertools
import gzip



from .biopython import FastqGeneralIterator

sys.dont_write_bytecode = True


#### USEFUL STUFF ####

def get_ambig_dna_values():
    '''
    Dict of IUPAC alphabets ambigious base calls
    '''
    ambiguous_dna_values = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "M": "AC",
        "R": "AG",
        "W": "AT",
        "S": "CG",
        "Y": "CT",
        "K": "GT",
        "V": "ACG",
        "H": "ACT",
        "D": "AGT",
        "B": "CGT",
        "X": "GATC",
        "N": "GATC",
        }
    return ambiguous_dna_values


def get_ambig_complement():
    '''
    Dict of IUPAC base call complement
    '''
    ambiguous_dna_complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "M": "K",
        "R": "Y",
        "W": "W",
        "S": "S",
        "Y": "R",
        "K": "M",
        "V": "B",
        "H": "D",
        "D": "H",
        "B": "V",
        "X": "X",
        "N": "N",
        }
    return ambiguous_dna_complement


# Line counter #

def file_len(fname):
    '''
taken from
http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python?lq=1
    '''
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# This just converts amig bases to bracketed bases. #
# Useful for making a 'grep'-able string #
def ambig_base_remover(in_seq):
    nuc_list = ['A', 'C', 'G', 'T']
    ambig_base_dict = get_ambig_dna_values()
    fixed_seq = ''
    for nuc in in_seq:
        if nuc in nuc_list:
            fixed_seq = fixed_seq + nuc
        if nuc not in nuc_list:
            actual_base = '[' + ambig_base_dict[nuc] + ']'
            fixed_seq = fixed_seq + actual_base
    return fixed_seq



# Simple reverse complementer #
def reverse_complement(seq):
    comp_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reversed_seq = seq[::-1]
    rev_comp_seq = ''
    for nuc in reversed_seq:
        rev_comp_seq = rev_comp_seq + comp_dict[nuc]
    return rev_comp_seq


# This generates  the files for appending #
def outfile_maker(files_path, sample_dict, dir_opt=None):
    if dir_opt is None:
        cur_dir = files_path
        fq_end = '.fastq'
    if dir_opt is not None:
        cur_dir = os.path.join(files_path, dir_opt)
        fq_end = '_' + dir_opt + '.fastq'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)
    for key in sample_dict.keys():
        cur_handle = os.path.join(cur_dir, (key + fq_end))
        with open(cur_handle, 'w+') as f:
            pass
    return cur_dir

# writes files and leaves open #
def outfile_opener(files_path, sample_dict, file_mode='w+', dir_opt=None):
    file_dict = {}
    if dir_opt is None:
        cur_dir = files_path
        fq_end = '.fastq'
    if dir_opt is not None:
        cur_dir = os.path.join(files_path, dir_opt)
        fq_end = '_' + dir_opt + '.fastq'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)
    for key in sample_dict.keys():
        cur_handle = os.path.join(cur_dir, (key + fq_end))
        if file_mode == 'w+':
            file_dict[key] = open(cur_handle, file_mode)
        if file_mode == 'wb':
            file_dict[key] = gzip.open((cur_handle + '.gz'), file_mode)
    return file_dict


# This makes the bc dicts #
def bc_dict_maker( bc_file, rev_comp_bcs):
    bc_dict = {}
    bc_len = []
    for line in open(bc_file, 'rU'):
        comb_id, bc = line.split()
        bc_len.append(len(bc))
        rev_bc = reverse_complement(bc)
        if rev_comp_bcs is True:
            bc_dict[comb_id] = rev_bc
        if rev_comp_bcs is False:
            bc_dict[comb_id] = bc
    return (bc_dict, max(bc_len))

# This just trims the read #
def trimmed_read(seq, qual, start_pos, end_pos, fp_len, bp_len, infile_extension):
    cur_seq = seq\
            [(start_pos + fp_len):(len(seq) - (end_pos + bp_len))]
    if infile_extension == 'fasta':
        cur_qual = None
        return (cur_seq, cur_qual)
    if infile_extension == 'fastq':
        cur_qual = qual\
            [(start_pos + fp_len):(len(seq) - (end_pos + bp_len))]
        return (cur_seq, cur_qual)


# Creates a dict with sampleID as key and
# list as item
def sample_index_maker(dict1, dict2):
    index_dict = {}
    sample_list = []
    for r1 in dict1.keys():
        for r2 in dict2.keys():
            sample_list.append(r1 + r2)
    for i in sample_list:
        index_dict[i] = []
        index_dict['Unassigned'] = []
    return index_dict

# this is what writes each of ind files by appending #
def ind_fastq_writer(cur_index, files_path, dir_opt=None):
    if dir_opt is None:
        cur_dir = files_path
        fq_end = '.fastq'
    if dir_opt is not None:
        cur_dir = files_path #+ '/' + dir_opt
        fq_end = '_' + dir_opt + '.fastq'
    for key, seq_list in cur_index.items():
        for seq_rec in seq_list:
            hd, seq, qual = seq_rec
            cur_handle = os.path.join(cur_dir, (key + fq_end))
            with open(cur_handle, 'a+') as f:
                f.write('@' + hd + '\n'+ seq + '\n+\n' + qual + '\n')

def ind_fastq_chunk_writer(cur_index, files_path, dir_opt=None):
    if dir_opt is None:
        cur_dir = files_path
        fq_end = '.fastq'
    if dir_opt is not None:
        cur_dir = files_path
        fq_end = '_' + dir_opt + '.fastq'
    for key, seq_chunk in cur_index.items():
        cur_handle = os.path.join(cur_dir, (key + fq_end))
        with open(cur_handle, 'a+') as f:
            f.write(str(seq_chunk))



def single_fasta_fastq(seq_rec, file_extension, out_file, out_qual_file, min_seq_length):

    header, cur_seq, cur_qual = seq_rec

    if file_extension == 'fasta':
        if len(cur_seq) < min_seq_length:
            return
        else:
            out_file.write('>' + header + '\n' + cur_seq + '\n')
            if out_qual_file is not None:
                out_qual_file.write('>' + header + '\n' + cur_qual + '\n')


    if file_extension == 'fastq':
        if len(cur_seq) < min_seq_length:
            return
        else:
            out_file.write('@' + header + '\n' + cur_seq \
            + '\n' + '+' +'\n' + cur_qual + '\n')




# Determines infile format #
def infile_extension(infile):
   with open(infile, 'rU') as f:
       first_char = f.read(1)
       if first_char == '@':
           return(str('fastq'))
       if first_char == '>':
           return(str('fasta'))

 # This creates a dict w/ with the seqs and a set w/ #
 # header info #
def read_set_maker( seq_batch):
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

def demult_header( r1_id, r2_id, seq_number, orig_header, keep_original_headers):
   if 'Unk' in (r1_id + r2_id):
       sample_id = 'Unassigned'
   if 'Unk' not in (r1_id + r2_id):
       sample_id = r1_id + r2_id
   if keep_original_headers is True:
       header = orig_header + ' ' + sample_id + '_' + str(seq_number)
   if keep_original_headers is False:
       header = sample_id + '_' + str(seq_number)
   return(sample_id, header)

 # This function takes the bc containing bits #
 # and matches them to the dicts #
def de_bc_loop( sub_seq, bc_dict):
    for key, bc in bc_dict.items():
        if bc in sub_seq:
            comb_id = key
            start_pos = len(bc)
            break
        if bc not in sub_seq:
            comb_id = 'Unk'
            start_pos = 0
    return (comb_id, start_pos)




# Batch Generator Function taken from Biopython website #
# http://biopython.org/wiki/Split_large_file #

def batch_iterator( iterator, batch_size) :
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = next(iterator)
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def mc_spinner(proc='\n'):
    sys.stdout.write(proc + '\n')
    sys.stdout.flush()
    spinner = itertools.cycle(['-', '/', '|', '\\'])
    while True:
        sys.stdout.write(next(spinner))
        sys.stdout.flush()
        time.sleep(0.1)
        sys.stdout.write('\b')
