#!/usr/bin/env python

"""
ParseSumaClust.py
Take from the output of sumaclust
and produce a fasta with the ;size=XXX field
@author: amnon
"""

__version__ = "0.9"

import argparse

from os import listdir
from os.path import isfile,join
import sys
import numpy as np
from cogent.parse.fasta import MinimalFastaParser
import re

def ParseSumaClust(inputfilename,outfilename):
	allfreqs={}
	allseqs={}

	# load the sequences
	infile=open(inputfilename)
	numreadsp=re.compile('count=(\d+)')
	clusteridp=re.compile('cluster=(\d+)')
	clusterscorep=re.compile('cluster_score=(\d+\.+\d+)')
	for seqid,seq in MinimalFastaParser(infile):
		newseqid=seqid
		res=numreadsp.search(newseqid)
		numreads=int(res.group(1))

		res=clusteridp.search(newseqid)
		clusterid=res.group(1)

		res=clusterscorep.search(newseqid)
		clusterscore=float(res.group(1))

		if allfreqs.has_key(clusterid):
			allfreqs[clusterid]+=numreads
		else:
			allfreqs[clusterid]=numreads

		# if it is the cluster center store the sequence
		if clusterscore==1:
			allseqs[clusterid]=seq
	infile.close()

	# now save the new fasta file
	numclusts=0
	with open(outfilename,'w') as outfile:
		for (clusterid,reads) in allfreqs.items():
			if allseqs.has_key(clusterid):
				outfile.write('>'+clusterid+';size='+str(reads)+'\n')
				outfile.write(allseqs[clusterid]+'\n')
				numclusts+=1
			else:
				raise NameError("Cluster with no center: "+clusterid)

	print("total clusters: "+str(numclusts))


def main(argv):
    parser=argparse.ArgumentParser(description='add size annotations to uparse fasta '+__version__)
    parser.add_argument('-i','--input',help='sumaclust output fasta')
    parser.add_argument('-o','--output',help='output fasta file name')
    args=parser.parse_args(argv)
    ParseSumaClust(args.input,args.output)
    
if __name__ == "__main__":
    main(sys.argv[1:])                
