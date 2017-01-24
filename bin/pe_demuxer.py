#!/usr/bin/python

from __future__ import division

from Mr_Demuxy import pe_demuxer_dist

def main():
    pe_class = pe_demuxer_dist.PEDemux()
    pe_class.main_loop()
if __name__ == '__main__':
    main()

