#!/usr/bin/python


from mr_demuxy import merged_demuxer_dist

def main():
    merged_class = merged_demuxer_dist.MergedDemuxer()
    merged_class.main_loop()

if __name__ == '__main__':
    main()
