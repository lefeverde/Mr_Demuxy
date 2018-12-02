#!/usr/bin/python


from mr_demuxy import pe_demuxer_dist
#from Mr_Demuxy import pe_demuxer_dist


def main():
    pe_class = pe_demuxer_dist.PEDemux()
    pe_class.main_loop()

if __name__ == '__main__':
    main()
