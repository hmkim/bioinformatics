import os, sys

from back_end import Downloader

from os.path import basename

# now you can call it directly with basename
import ntpath


#tax_s=[6248  ,51022 ,51031 ,6252  ,102285,6216  ,85433 ,51028 ,6206  ,6204  ,60517 ,31281 ,58839,] 
tax_s=[6248]

tree_map={}

from Bio import SeqIO
from Bio import AlignIO

import subprocess

import argparse

from Bio.Align.Applications import ClustalOmegaCommandline 

def select_template(fa_file):
    max_cluster_weight = 0
    template = None
    template_seq = None
    for index, record in enumerate(SeqIO.parse(fa_file, "fasta")):
        record_des_s = (record.description).split(' ')
        if 'cluster_center=True;' in record_des_s:
            id_ = record_des_s.pop(0)
            dict_ = getValue(record_des_s)

            if not template:
                max_cluster_weight = dict_['cluster_weight']
                template_seq = record.seq
                template = id_
            if max_cluster_weight < dict_['cluster_weight']:
                max_cluster_weight = dict_['cluster_weight']
                template_seq = record.seq
                template = id_
    #print (max_cluster_weight, template)
    return (max_cluster_weight, template, template_seq) 


def find_oligo(seqname, seq):
    seq_len = len(seq)
    print (seqname, seq)
    sys.exit()
    k = seq_len - 100

    flag = 0

    oligo_length = 20

    while (k > oligo_length + 50):

        candidate = seq[k - oligo_length:k]
        print(candidate)
#        # If the criteria is not met, try again in the next 'oligo_length' region
        if (flag < 5):
            k = k - 1
        else:
            break

def getValue(list_):
    result = dict()
    for line in list_:
        key = line.split('=')[0]
        value = line.split('=')[1]
        value = value.replace(';','')
        result[key] = value
    return result 

def tree_draw_newick(aligned):
    alignment_overlapped= AlignIO.read(aligned, "fasta")
    
def align(in_file,out_file):
    if not in_file or not out_file:
        print ("ERROR! no file for align")
        sys.exit()
    clustal_commandline = {'infile': in_file,
                           'outfile': out_file,
                           'outfmt': 'fasta',
#                           'iterations': 3,
                            'cmd': '/Users/brandon/Downloads/bio/preprocess/tools/clustalo',
                           'force': True}

    try:
        clustalomega_cline = ClustalOmegaCommandline(**clustal_commandline)
    except Exception as e:  # This will catch all the major errors
        print ('Error in command directed to ClustalO. Check the command entered!')
        print ('Error produed:\n%s' % e)
        raise

    try:
        clustalomega_cline()
    except Exception as e:  # This will catch all the major errors
        print ('Error when running Clustal Omega. Make sure Clustal Omega is installed and working properly')
        print ('Error produed:\n%s' % e)
        raise

def gb_to_fasta(gbfile,fastafile):
    fw = open(fastafile, 'w')
    for index, record in enumerate(SeqIO.parse(gbfile, "genbank")):

        cons_seq = ''
        cons_seq_string = []
        for i,f in enumerate(record.features):
            if (i == 0): ## pass 'source' type
                continue

            strand = str(f.strand)

            product = f.qualifiers.get("product", ["-"])[0]
            note = f.qualifiers.get("note", ["-"])[0]

            try:
                seq = f.location.extract(record).seq
            except ValueError:
                seq = record.seq[f.location.start.position:f.location.end.position]
                print("%s:ValueError"%record.id)

            # check strand
            if strand == "-1":
                seq = seq.reverse_complement()
            else:
                pass

            #if '<' in str(f.location) or '>' in str(f.location):
            #    continue
            
            #print ("%s\t%s\t%s\t%s\t%s\t%s"%(record.id, f.type, product, note, f.location,  seq)) 
            cons_seq_string.append('product:' + product + ' ' + str(f.location))
            cons_seq += seq

        if cons_seq:
            print (">%s %s\n%s"%(record.id," ".join(cons_seq_string), cons_seq))    
            fw.write (  (">%s %s\n%s\n"%(record.id," ".join(cons_seq_string), cons_seq))    )
    fw.close()
            



def process(taxid):

    search_term = 'txid' + str(taxid) + '[Organism] AND ("ITS1"[All Fields] OR "ITS 1"[All Fields] OR "internal transcribed spacer 1"[All Fields]) AND ("ITS2"[All Fields] OR "ITS 2"[All Fields] OR "internal transcribed spacer 2"[All Fields]) NOT mitochondrial[All Fields]'
    print (search_term)
    sys.exit()

    gbfile = 'data/txid' + str(taxid) + '.gb'

    # download
#    if not (os.path.isfile(gbfile)) or os.stat(gbfile).st_size == 0:
#        dler = Downloader('nucleotide', search_term, gbfile, 0, email="mail@gmail.com")
#        dler.run_everything()

def msf_to_fasta_using_perl(msf,fasta):
    cmd = 'source /Users/brandon/perl5/perlbrew/etc/bashrc; perl msf_to_fa.pl %s %s'%(msf, fasta)
    print (cmd)
    subprocess.call(cmd, shell=True)
   
def vsearch_for_fasta(fasta,output):
    cmd = 'vsearch --threads 4 --cluster_fast %s --id 0.70 --clusters %s'%(fasta,output)
    print (cmd)
    subprocess.call(cmd, shell=True)

def main():

#    for tax in tax_s:
#        process(tax)

    parser = argparse.ArgumentParser(description='Oligo Design Runner')
    parser.add_argument('-f', '--fasta', help='fasta file')
    parser.add_argument('-g', '--genbank', help='genbank file')
    parser.add_argument('-a', '--aligned', help='aligned file')
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    args = parser.parse_args()

    gb_file = args.genbank
    fa_file = args.fasta 

    if not args.fasta:
        print ("no fasta")
        sys.exit()
    else:
        fa_file = args.fasta

    if not args.aligned:
        basepath = os.path.dirname(fa_file)
        base = os.path.basename(fa_file)
        basename = os.path.splitext(base)[0]
        aligned_file = basepath + basename + '.msf'        
    else:
        aligned_file = args.aligned
  
    # collect sequence
      
    # filtering using vsearch, sumaclust
      
    # align sequence
    #align(fa_file,aligned_file)

    
    #results=tree_draw_newick(fa_file)
    
    #msf_to_fasta_using_perl(aligned_file,fa_file)
    #gb_to_fasta(gb_file,fa_file)

    #vsearch_file = basepath + basename + '.vsearch'
    #vsearch_for_fasta(fa_file, vsearch_file)

    #print (vsearch_file)

    # select record in cluster
    (weight,temp_id,temp_seq) = select_template(fa_file)

    # find oligo candidate in template sequence 
    find_oligo(temp_id, temp_seq)
        
   

if __name__ == '__main__':
    main()    

