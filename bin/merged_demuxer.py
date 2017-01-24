#!/usr/bin/python

from __future__ import division

from Mr_Demuxy import merged_demuxer_dist

def main():
    merged_class = merged_demuxer_dist.MergedDemuxer()
    merged_class.main_loop()
    
if __name__ == '__main__':
    main()

