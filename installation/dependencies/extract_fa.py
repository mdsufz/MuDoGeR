#!/usr/bin/env python 

# this script takes three inputs files.
#                       1) list of ids
#                       2) multi-fasta file
#                       3) output file with selected fasta

import os
import sys
from itertools import groupby
def get_ids(id_file):
        # get ids from id_file and put into dictionary
        # return dic
        l="a"
        f=open(id_file,"r")
        ids = f.readlines()
        dic = {}
        for id in ids:
                dic[id.replace(">","").split(" ")[0].strip()]=""

        f.close()

	#print dic
        return dic

def extract_fa(id_dic,fa_in,fa_out,type):
        # given the id dic and the fasta file, extracts the desired fasta seqs
        pass
        l="a"
        #f=open(fa_in,"r")
        fout = open(fa_out,"w")
        dic = {}
	#print "type",type
	
	#print id_dic
	if type == 1:	
        	f=open(fa_in) # opens fasta file
        # ditch the boolean (x[0]) and just keep the header or sequence since
         # we know they alternate.
        	f_iterator = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
        	for header in f_iterator:
               		header = header.next()[1:].strip()                   # drop the ">"   # <-- optional
                	seq = "".join(s.strip() for s in f_iterator.next())  # join all sequence lines to one.
                	seq = str(seq)                                       # cast sequence to string  #  
                	#header = str(header)                                 # cast header to string    #  for troubleshooting purpose
                	aux = header.replace(">","").split(" ")[0]
		#	aux = aux.replace("_",".")
		#	print aux
			
			if aux in id_dic:
        			#print "foi"               		
				fout.write(">"+header+"\n")
                        	fout.write(seq+"\n")

	        f.close()

	else:
		f=open(fa_in,"r")
		l="a"
		print len(id_dic)
		while True:
			l=f.readline()
			if not l:break
			else:
				if ">" in l:
					#print l
					header = l.strip()
					#print header
					seq = f.readline().strip()
					aux=header.replace(">","")#.split(" ")[0]
				#	aux=aux.replace("_",".")
				#	print aux	
					#print aux
					
					 
					aux=header.split()[0].replace(">","")
					#print aux
					if aux in id_dic:
						#print "foi"
						fout.write(header+"\n")
                                               	fout.write(seq+"\n")
					 
					#aux=header.split()[0].replace(">","")
					#print aux
				################ sweeping thorugh the dic #####
					#for key,value in id_dic.iteritems():
						#print header
					#	if key in header:
					#	key=key[1:]
						#print key, aux
						#dsfdas= raw_input()
					#	if key == aux:
					#		print "foi"
					#		fout.write(header+"\n")
					#		fout.write(seq+"\n")		
					#
				############################################
def check_type(fa_in):
	c = 0
	list = []
	f=open(fa_in,"r")
	l="a"
        while True:
        	l=f.readline()
                if not l:break
                else:
                	if ">" in l:
				list.append(l)
				c = c +1 
				if c == 5: break
	if ">" in list[0]:
		if ">" in list[2]:
			f.close()
			return 2	# if first and third lines have headers, this is a long line multi-fasta file
		else:
			f.close()
			return 1	# else, its a "broken" line fasta file
	
	

def main(id_file,fa_in,fa_out):
        
	#type = check_type(fa_in)	# check what type of multi fasta file you have: long line or broken line
	type =1	
	# extract dictionary
	id_dic = get_ids(id_file)
	#print id_dic
        # dumps to the outputfile
        extract_fa(id_dic,fa_in,fa_out,type)


# get inputs from user keyboard
id_file,fa_in,fa_out = sys.argv[1], sys.argv[2], sys.argv[3]
main(id_file,fa_in,fa_out)
