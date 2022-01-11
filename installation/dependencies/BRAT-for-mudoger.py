import os
import math
import subprocess
import sys
import commands
import os.path


for i in range(1):	# HEADER
	# INPUTS: bins_folder, libraries_folder, output_folder

	#algorithm:
	#1) check if all bins in bins_folder are indexed with bowtie2-build
		# 1.1) if not, run bowtie2-build
		
	#2) run bowtie2 between all the bins and all the libraries, dumping the sam file to output_folder

	#3) run samtools view -F 4 in all sam files in output_folder and dump to a text file in output_folder

	#4) process text file in output_folder and format it to a tsv table -> brat

	a=1

def module_load():
	os.system("module load bowtie2 samtools")

def check_bins(bins_folder): # step 1
	bins_a=os.listdir(bins_folder)
	bins=[]
	for b in bins_a:
		bins.append(bins_folder+"/"+b)
	for b in bins:
		if os.path.isfile(b.replace(".fa",".1.bt2")):
			pass
		else:
			status,aux=commands.getstatusoutput("bowtie2-build -q "+b+" "+b.replace(".fa",""))
	print "-----> bins are okay"
	

	
def run_mapping(bins_folder, libraries_folder, output_folder): # step 2
	threads=str(2)
	# get bins
	bins_a=os.listdir(bins_folder)
	bins=[]
	for b in bins_a:
		if "fa" in b:
			bins.append(bins_folder+"/"+b.replace(".fa",""))
	# get libs
	libs_a=os.listdir(libraries_folder)
	libs=[]
	for l in libs_a:
		if "samp"  in l:
			libs.append(libraries_folder+"/"+l)
	
	#map
	for b in bins:
		for l in libs:
			sam=output_folder+"/"+b.split("/")[-1]+"__"+l.split("/")[-1].replace(".fa","")+".sam"
			#print b,l,sam
			status,aux=commands.getstatusoutput("bowtie2 --threads "+threads+" -x "+b+" -f "+l+" -S "+sam)
			#print status, aux
	print "-----> mapping was done"

def run_samtools(output_folder):
	mappings_a=os.listdir(output_folder)
	mapping=[]
	map_dic={}
	for m in mapping_a:
		mapping.append(output_folder+"/"+m)
		status,aux=commands.getstatusoutput("samtools view -F 4 "+output_folder+"/"+m+" ")
				
		
		
		
			
def main (bins_folder, libraries_folder, output_folder):
	module_load()
	#1
	check_bins(bins_folder)
	#2
	run_mapping(bins_folder,libraries_folder,output_folder)
	#3
	run_samtools(output_folder)
	
	
	
	
	
bins_folder="/data/msb/tools/scripts/BRAT/test/bins"
libraries_folder="/data/msb/tools/scripts/BRAT/test/libs"
output_folder="/data/msb/tools/scripts/BRAT/test/output"
#bins_folder, libraries_folder, output_folder = sys.argv[1],sys.argv[2],sys.argv[3]
main(bins_folder, libraries_folder, output_folder)
