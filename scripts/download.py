import os

from back_end import Downloader

#tax_s=[51022, 6248  ,51022 ,51031 ,6252  ,102285,6216  ,85433 ,51028 ,6206  ,6204  ,60517 ,31281 ,58839] 
tax_s=[51028]


def process(taxid):

    #search_term = 'txid' + str(taxid) + '[Organism] AND ("ITS1"[All Fields] OR "ITS 1"[All Fields] OR "internal transcribed spacer 1"[All Fields]) AND ("ITS2"[All Fields] OR "ITS 2"[All Fields] OR "internal transcribed spacer 2"[All Fields]) NOT mitochondrial[All Fields]'
    search_term = 'txid' + str(taxid) + '[Organism] AND (internal transcribed spacer)'

    print ("esearch -db nuccore -query '%s' | efetch -format fasta"%search_term)

    
    gbfile = 'data/txid' + str(taxid) + '.gb'

    # download
#    if not (os.path.isfile(gbfile)) or os.stat(gbfile).st_size == 0:
#        dler = Downloader('nucleotide', search_term, gbfile, 0, email="brandon.kim.hyunmin@gmail.com")
#        dler.run_everything()


def main():

    for tax in tax_s:
        process(tax)


if __name__ == '__main__':
    main()    

