#https://github.com/brantfaircloth/sliding_primer/blob/master/sliding.py
#https://github.com/savill88/PriMSA/blob/b1a2adc00f31fa37f6212d8205289754b8b5e3c5/django_PriMSA/primsa/views.py
#https://github.com/semal/GenoPrime2/blob/f0e7fcb4832072440b2e460f58bd01fbe73b8ab5/genoprime/primer.py
#https://libnano.github.io/primer3-py/quickstart.html#primer-design
#http://public.lanl.gov/jgans/tntblast/tntblast_doc.html
#https://github.com/allista/DegenPrimer
#https://github.com/joshquick/primal
#https://github.com/qPCR4vir/VisualOliDeg

import epride as ep

from Bio.Align import MultipleSeqAlignment


from Bio.Align.Applications import ClustalOmegaCommandline 

from Bio import SeqIO, AlignIO

def align(in_file,out_file):
    clustal_commandline = {'infile': in_file,
                           'outfile': out_file,
                           'outfmt': 'fasta',
#                           'iterations': 3,
                            'cmd': '/Users/brandon/Downloads/dev_AutoMSA/tools/clustalo',
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

    
def process_msf(aligned):
    #alignment = AlignIO.read(align, "clustal")
    alignment_overlapped= AlignIO.read(aligned, "fasta")
    print (alignment_overlapped)

def process(fasta):
    # Cluster the sequences in test_align.fasta
    test_clusters = ep.Primers(fasta, id_val=0.8, primer_length=20)

    test_clusters.coverage_info(id_penalty=1, primer_penalty=1)
    
    for i in test_clusters[0]:
        print (i)

def design(record,  start = 0, stop = 150):
    while stop < len(record.seq):
        seq_slice = str(record.seq[start:stop])
        print (seq_slice)
        stop = start + 150



def main():
    #output = open('primersDesigned.txt','w')
    #count = 0

    fasta_file = '/Users/brandon/Downloads/dev_AutoMSA/data/txid6248.gb.fa'
    #in_msf = '/Users/brandon/Downloads/dev_AutoMSA/data/tax_6216.gb.97.msf'
    aligned_file = '/Users/brandon/Downloads/dev_AutoMSA/data/txid6248.gb.aligned.fa'

    align(fasta_file,aligned_file)

    #process(in_fasta)
    process_msf(aligned_file)

    # read the sequence
#    for record in SeqIO.parse(open(in_fasta, 'r'), 'fasta'):
#        design(record)
#        sys.exit()
#    output.close()

if __name__ == '__main__':
    main()

