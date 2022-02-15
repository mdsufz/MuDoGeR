#!/usr/bin/env python 

# this script takes two inputs files.
#                       1) multi fasta file 
#                       2) location to dump the splitted fasta files (one file for each sequence)

import os
import sys
from itertools import groupby

def split_fa(fa_in,out_base):
	l="a"
	c=0
	f=open(fa_in,"r")	
	l="a"
	while True:
		l=f.readline()
		if not l:break
		else:
			if ">" in l:
				header = l.strip()
				seq = f.readline().strip()
				fout = open(out_base+"-"+str(c)+".fa","w")
				fout.write(header+"\n")
				fout.write(seq+"\n")
				fout.close()
				c=c+1
				 

def main(fa_in,out_base):
        
        
        split_fa(fa_in,out_base)


# get inputs from user keyboard
fa_in,out_base = sys.argv[1], sys.argv[2] 
main(fa_in,out_base)
