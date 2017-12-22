#from Bio import Entrez

from back_end import Downloader

import sys

def main():
    search_term = sys.argv[1]
    outfile = sys.argv[2]


    dler = Downloader('nucleotide', search_term, outfile, 0)
    dler.run_everything()

if __name__ == '__main__':
    main()
