import os
import sys
import commands
import random



def dict_replenish(list,dic_features):
	dict_bins={}
	for item in list:
		dict_bins[item]=dic_features[item]

	return dict_bins

def select_bin(bins_quality_dic):
	#print '######################### SELECT BINS'
	dict_bins=bins_quality_dic
	#for bin in bins:
	#	dict_bins[bin]=bins_quality[bin]
	#	print bin,bins_quality[bin]

	# by completeness
	comp=0
	best=[]
	#aux_dic = dict_bins
	for k,v in dict_bins.iteritems():
		if v[0]>comp:
			comp=v[0]
			best=[]
			best.append(k)
		elif v[0] == comp:
			best.append(k)
			pass
	if len(best) > 1:
		# by contamination
		dict_bins = dict_replenish(best,dict_bins)  # restore dictionary
		best = []
		cont=100
		for k,v in dict_bins.iteritems():
			if v[1]<cont:
				cont=v[1]
				best=[]
				best.append(k)
			elif v[1] == cont:
				best.append(k)
				pass
		if len(best) > 1:

			# by strain heterogeneity
			dict_bins = dict_replenish(best,dict_bins)  # restore dictionary
			best = []
			strhet=100
			for k,v in dict_bins.iteritems():
				if v[2]<strhet:
					strhet=v[2]
					best=[]
					best.append(k)
				elif v[2] == strhet:
					best.append(k)
					pass
			if len(best) > 1:

				# by n50
				dict_bins = dict_replenish(best,dict_bins)  # restore dictionary
				best = []
				n50=0
				for k,v in dict_bins.iteritems():
					if v[3]>n50:
						n50=v[3]
						best=[]
						best.append(k)
					elif v[3] == n50:
						best.append(k)
						pass
				if len(best) > 1:

					# by number of contigs
					dict_bins = dict_replenish(best,dict_bins)  # restore dictionary
					best = []
					numcontigs=100000
					for k,v in dict_bins.iteritems():
						if v[4]<numcontigs:
							numcontigs=v[4]
							best=[]
							best.append(k)
						elif v[4] == numcontigs:
							best.append(k)
							pass
					
					if len(best) > 1:

						# by number of nucleotides
						dict_bins = dict_replenish(best,dict_bins)  # restore dictionary
						best = []
						size_genome=1
						for k,v in dict_bins.iteritems():
							if v[5]> size_genome:
								size_genome=v[5]
								best=[]
								best.append(k)
							elif v[5] == size_genome:
								best.append(k)
								pass
	if len(best) > 1:							 # if no bin is chosen based on the criteria above, pick randomly
		#for k,v in dict_bins.iteritems():
		#	print k,v
		representative = random.choice(best)
		#print 'representative:',representative
	else: 
		representative = best[0]

	return representative
	
sep="\t"
#print "CLUSTER"+sep+"SUB-CLUSTER"+sep+"BIN"
cluster_folder=sys.argv[1]
checkm_file=sys.argv[2]
n50_file=sys.argv[3]
numcontigs_size_file=sys.argv[4]
cluname=cluster_folder.split()[-1].split('-')[-1]
def remove_duplicate(duplicate):
	final_list = []
	for num in duplicate:
		if num not in final_list:
			final_list.append(num)
	return final_list
taxonomy=cluster_folder+"/taxonomy"
t=open(taxonomy,"r")
tax=t.readline()
t.close()


# get'em clusters
#summary=cluster_folder+"/ani_split_95/cluster_summary.tsv"
summary=cluster_folder+"/ani-split-95-ctonly/cluster_summary.tsv"
f=open(summary,"r")
header=f.readline()
clu_dic={}
while True:
    l=f.readline()
    if not l:   break
    else:
       # print (l)
      #  print (l.split())
        clu=l.split()[0]
        if clu not in clu_dic:
            clu_dic[clu]=[l.split()[1]]
            #clu_dic[clu].append(l.split()[2])
        else:
            if l.split()[1] not in clu_dic[clu]:
                clu_dic[clu].append(l.split()[1])
            if l.split()[2].strip() not in clu_dic[clu]:
                clu_dic[clu].append(l.split()[2].strip())
f.close()

#for k,v in sorted(clu_dic.iteritems()):
#	print k,v

#print '###### SUB CLUSTERS EXTRACTED - NOW LETS SELECT REPRESENTATIVES OF EACH SUBCLUSTER'

cluster_dic_bins={}
for k,v in clu_dic.items():
	#list_bins_clu=[k]
	cluster_dic_bins[k]=[]
#	print 'K:',k
	
	for bin in v:
		# get completeness, contamination and strain heterogeneity
		status,output=commands.getstatusoutput("grep -w "+bin.replace(".fa","")+" "+checkm_file+" | cut -f2,3,4")
		comp,cont,strhet=output.split("\t")
		
		# get n50
		#print bin
		status,n50=commands.getstatusoutput("grep -w "+bin.replace(".fa","")+" "+n50_file+" | cut -f2 ")
		
		# get number of contigs and size of bin
		status,output=commands.getstatusoutput("grep -w "+bin.replace(".fa","")+" "+numcontigs_size_file+" | cut -f2,3")
		num_contigs,size=output.split("\t")
		

		#print ("cluster-"+cluname+"\tsubcluster-"+k+"\t"+bin) #+"\t"+tax.strip())

		#print 'size',size
		#list_bins_clu.append([bin,float(comp),float(cont),float(strhet),int(n50),int(num_contigs),int(size)])
		cluster_dic_bins[k].append([bin,float(comp),float(cont),float(strhet),int(n50),int(num_contigs),int(size)])

		

    	#list_bins_clu.sort(key=lambda x:(x[1], -x[2], x[3],), reverse=True) # sort list of bins by descending completeness and ascending contamination
    	
		#print len(list_bins_clu)	

for k,v in cluster_dic_bins.iteritems():
	#print k
	subcluster_dic={}
	for bin in v:
	#	print 'bin:',bin
	#	print 'v:',v[0]
		subcluster_dic[bin[0]]=[bin[1],bin[2],bin[3],bin[4],bin[5],bin[6]]
	for t,b in subcluster_dic.iteritems():
		#print t,b
		pass
	rep= select_bin(subcluster_dic)
	print cluname+sep+k+sep+rep.replace(".fa","")
sys.exit()
#print list_bins_clu

#	list_bins_clu = remove_duplicate(list_bins_clu) 
#	#print len(list_bins_clu)
#	if len(list_bins_clu) > 1:
#		for bin in list_bins_clu:
#			if list_bins_clu[0][1]==list_bins_clu[1][1] and  list_bins_clu[0][2] ==  list_bins_clu[1][2] and list_bins_clu[0][3] ==  list_bins_clu[1][3] and list_bins_clu[0][4] ==  list_bins_clu[1][4] and  list_bins_clu[0][5] ==  list_bins_clu[1][5]:
#				print "we have problems here >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#				for bin in list_bins_clu:
#					print bin
#   	else:
#		print "representative:",list_bins_clu[0]
	




print "--------------------------------------------------------------------------------------------------------------------------------------------------"
