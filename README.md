


Mr\_Demuxy 
================



A Demultixper for Combinatorially Multiplexed Next Generation Reads

Author: Daniel E. Lefever

Author-email
[lefeverde<span>@</span>gmail<span>.</span>com](mailto:lefeverde%40gmail.com)

License: MIT


Introduction:
-------------

Mr. Demuxy is a fairly simple, self-contained Python program that makes
demultiplexing variable length, combinatorially barcoded NGS reads as
painless as possible. See below for more details. This program was
written with the combinatorial barcoding strategy developed by Travis
Glenn in mind \[Adapterama Ref\], however, it would work on any type of
data multiplexed with combinatorial barcodes. In short, this approach
uses two pairs of barcodes denoted as the internal and external barcodes
(also referred to as tags and indexes). The idea is to pool many
different projects onto one MiSeq/HiSeq run cheaply by having external
Illumina barcodes to differentiate the pooled projects and internal
‘homebrewed’ variable length barcodes to distinguish samples within each
project. The internal barcoding strategy is based on the layout of a
96-well plate: the forward tag corresponds to row letter (A-H) and the
reverse tag corresponds to column number (1-12). This program uses both
tags to come up with the sample ID.



Dependencies:
-------------

-   Python 2.7



Installation:
-------------

    and unarchive the folder. To unarchive the folder, try clicking
    on it. If that doesn’t work, open terminal and change directory to

        $ tar -xvzf Mr_Demuxy-x.y.z-dist.tar.gz

-   **NOTE**: The ‘x.y.z’ above denotes version number. Change that to
    whatever version is downloaded. Also, [relevant
    xkcd](http://xkcd.com/1168/)

1.  Then open up terminal and change directory to wherever the
    unarchived folder lives. Then, install it to your local
    site-packages directory like so:

        $ python setup.py install --user

2.  Restart terminal. You should now be able to run the program from
    anywhere in your directory.


### Can I use pip install?

Yes, but for whatever reason, pip install does not execute some extra
stuff I included in the setup script. The extra stuff basically adds a
line to your .bash\_profile which puts the local script directory in the
path. This is so you can call the scripts from anywhere. Also, it will
remove some stuff leftover from previous Mr. Demuxy versions. In order
to do this manually, you’ll need to add this to your .bash\_profile:

    export PATH=/full_path_to_scripts/:$PATH




Usage:
------

Mr. Demuxy contains two main scripts; merged\_demuxer.py and
pe\_demuxer.py for merged or paired-end reads, respectively. Basic usage
for each script are outlined below.

**Demultiplex Merged Reads (merged\_demuxer.py)**

    usage:
    $ merged_demuxer.py -f forward_bcs.txt -r reverse_bcs.txt -i example_merged.fastq  -o demultiplexed_seqs.fasta

    This program demultiplexes merged combinatorially barcoded reads and can
    output the demultiplexed reads a number of different ways.

    required arguments:
      -f                    forward bcs file
      -r                    reverse bcs file
      -i                    input merged fastq file

    optional arguments:
      -o                    name of output file or folder (default:
                            merged_demuxer_out)
      --forward_primer_len
                            base number of forward primer/adapter (default: 0)
      --reverse_primer_len
                            base number of reverse primer/adapter (default: 0)
      --no_rev_comp_reverse_bcs
                            prevents rev bcs from reverse complementing
      -h, --help            show this help message and exit
      --keep_original_headers
                            pass this arg to keep original headers in
                            demultiplexed fasta/fastq file
      --write_qual          write a qual file containing quality scores
      --individual_fastq    demultiplex by writing a fastq file for each sample ID
      --min_seq_length      minimum length of seq to be retained (default: 20)
      --file_extension      format of the output file (default: fasta)

**Demultiplex Paired-End Reads (pe\_demuxer.py)**

    usage:
    $ pe_demuxer.py -r1 r1_seqs.fastq -r2 r2_seqs.fastq -r1_bc r1_bc_file.txt -r2_bc -r2_bc_file.txt

    This program splits combinatorially multiplexed paired-end seqs and outputs a
    folder for the R1 reads and another for the R2 reads. The two folders contain
    one fastq file per sample id.

    required arguments:
      -r1_bc                BC file for R1 reads
      -r2_bc                BC file for R2 reads
      -r1                   R1 fastq file
      -r2                   R2 fastq file

    optional arguments:
      -h, --help            show this help message and exit
      -o                    name of output file or dir
      --batch_length        amount of seqs to be put in a batch only used when
                            check_headers is passed (Default: 1500000)
      --forward_primer_len
                            base number of primer/adapter on R1 seq (default: 0)
      --reverse_primer_len
                            base number of primer/adapter on R2 seq (default: 0)
      --keep_original_headers
                            pass this arg to keep original headers in
                            demultiplexed fasta/fastq file
      --min_seq_length      minimum length of seq to be retained (default: 20)
      --check_headers       ensures the R1 and R2 read match before demultiplexing
                            WARNING: this is super slow



Usage Examples: merged\_demuxer.py
----------------------------------

**Demultiplex merged fastq reads with output being in fasta format:**

    $ merged_demuxer.py -i example_merged.fastq -f forward_bcs.txt -r reverse_bcs.txt -o example_demult.fasta

**Do the same as above and get a seperate qual file:**

    $ merged_demuxer.py -i example_merged.fastq -f forward_bcs.txt -r reverse_bcs.txt -o example_demult.fasta --write_qual

**Demultiplex and keep the original fastq headers:**

    $ merged_demuxer.py -i example_merged.fastq -f forward_bcs.txt -r reverse_bcs.txt -o example_demult.fasta --keep_original_headers

**Demultiplex merged reads and remove the 17 bp forward primer and 22 bp
reverse primer:**

    $ merged_demuxer.py -i example_merged.fastq -f forward_bcs.txt -r reverse_bcs.txt -o example_demult.fasta --forward_primer_len 17 --reverse_primer_len 22

**Do the same as above except have the output be in fastq format:**

    $ merged_demuxer.py -i example_merged.fastq -f forward_bcs.txt -r reverse_bcs.txt -o example_demult.fastq --forward_primer_len 17 --reverse_primer_len 22 --file_extension fastq

**Demultiplex reads and split into one fastq file per sample ID:**

    $ merged_demuxer.py -i example_merged.fastq -f forward_bcs.txt -r reverse_bcs.txt -o example_demult_ind_fastq --individual_fastq

**Get individual fastq files, chop off the 17 bp forward primer and 22
bp reverse primer, and keep reads that are only longer than 100 bp:**

    $ merged_demuxer.py -i example_merged.fastq -f forward_bcs.txt -r reverse_bcs.txt -o example_demult_ind_fastq --individual_fastq --forward_primer_len 17 --reverse_primer_len 22 --individual_fastq --min_seq_length 100



Usage Examples: pe\_demuxer.py
------------------------------

**Demultiplex paired-end fastq files and write individual fastqs for
each sample in their respective directories:**

    $ pe_demuxer.py -r1 r1_seqs.fastq -r2 r2_seqs.fastq -r1_bc r1_bc_file.txt -r2_bc -r2_bc_file.txt

**Do the same as above and remove a 10 bp primer on R1 reads and a 16 bp
primer on R2 reads:**

    $ pe_demuxer.py -r1 r1_seqs.fastq -r2 r2_seqs.fastq -r1_bc r1_bc_file.txt -r2_bc -r2_bc_file.txt --forward_primer_len 10 --reverse_primer_len 16

**Demultiplex paired-end reads in which the read order may be different
and/or the files contain un-paired reads:** **WARNING:** using this
option is extremely slow and not recommended. Illumina sequencers should
produce files with reads in the same order and discard un-paired reads.
I didn’t realize that until after added this feature so I decided to
include it anyways.

    $ pe_demuxer.py -r1 r1_seqs.fastq -r2 r2_seqs.fastq -r1_bc r1_bc_file.txt -r2_bc -r2_bc_file.txt --check_headers



Input Files Information:
------------------------

This section explains more about the input files.

**Forward or R1 barcode file** The forward barcode file, as you might
expect, is the file containing the forward barcodes and letter (or
whatever the barcode corresponds to). You will probably need to make
your own forward barcode file, and the formatting is important. The
easiest way to make one, is to open a new sheet in excel and copy/paste
the barcode and letter pairs with the first column being the letter and
the second containing the barcode. It should look like this:

    A   GGTAC
    B   CAACAC
    C   ATCGGTT
    D   TCGGTCAA

-   **Note:** Both barcode files do not have any header information. I’m
    not sure if it will screw things up, but it’s safest to not
    include it. Once you’re done inputing the forward barcodes, save the
    excel sheet as a ‘tab delineated text’ file. Be sure to name it
    something so you know what it is.

**Reverse or R2 barcode file** The reverse barcode file is just like the
forward barcode file, except for the reverse barcodes. It should look
like this:

    1   AGGAA
    2   GAGTGG
    3   CCACGTC
    4   TTCTCAGC

-   **Note:** This program reverse complements (RC) the reverse barcodes
    when demultiplexing merged reads instead of RCing the whole
    sequence. This makes it go faster, but if your reverse barcodes do
    not need to be RC’d, pass this:

        --no_rev_comp_reverse_bcs



Output File Information:
------------------------

The output file/files contain reads that have been demultiplexed, and
optionally trimmed of primers. This program automatically removes the
barcode from each read and will remove the primers if length is
specified by the –forward\_primer\_len and –reverse\_primer\_len
arguments. The primers are not actually matched and are instead removed
by trimming a set number of bases after the barcode. This is much
quicker than actual matching.

For example, if the forward primer was 17 bp long, and the reverse 15,
pass:

    --forward_primer_len 17 --reverse_primer_len 15

If these options are not passed, it defaults to 0.

The merged\_demuxer.py script can also create individual fastq files
split by each sample. To do this, simply pass the –individual\_fastq
argument. The -o argument will become the name of the output directory
with one fastq per sample. The dir structure is as follows:

    output_dir
      |---<sample_id_1>.fastq
      |---<sample_id_2>.fastq
      .
      .
      |---<sample_id_n>.fastq

pe\_demuxer.py will only output individual fastq files. The R1 and R2
files are separated into two folders each containing individual fastq
files split by sample. This is the structure:

    output_dir
      |---R1
      |     |---<sample_id_1>_R1.fastq
      |     |---<sample_id_2>_R1.fastq
      |     .
      |     .
      |     |---<sample_id_n>_R1.fastq
      |
      |
      |---R2
            |---<sample_id_1>_R1.fastq
            |---<sample_id_2>_R1.fastq
            .
            .
            |---<sample_id_n>_R1.fastq

The merged\_demuxer.py script offers some control over how the headers
are formatted. If there’s some specific way that would be helpful,
please let me know and I’ll include it in future releases. Anyways, here
are the current options for the header format:

    HEADER OPTIONS:

    -Demultiplexed ID only-
    --FASTA Format--

    >A6_1
    AATATT....CACGCAGA

    --FASTQ Format--
    @A6_1
    AATATT....CACGCAGA
    +
    GGGGGG...@GGFGGFC<

    -Original header with demultiplexed ID-

    --FASTA Format--
    >M02849:...:1886 A6_1
    AATATT....CACGCAGA

    --FASTQ Format--
    @M02849:...:1886 A6_1
    AATATT....CACGCAGA
    +
    GGGGGG...@GGFGGFC<



Misc
----

This program was written and predominantly tested using a macbook pro
2011 with 8 GB memory and a 2.4 GHz Intel Core i7 processor running OS X
El Capitan using python 2.7.10. With this setup, both merged (written to
single or individual files) and paired-end demultiplexing of \~7,000,000
reads took around 3-4 minutes.



References:
-----------

*the fastq parser was taken from biopython*

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke,
A., … de Hoon, M. J. (2009). Bio\$ python: freely available \$ python
tools for computational molecular biology and bioinformatics.
Bioinformatics, 25(11), 1422-1423. doi:10.1093/bioinformatics/btp163

*Adapterama ref*
