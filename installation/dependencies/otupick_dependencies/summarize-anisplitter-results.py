#!/usr/bin/python3
# Version 1.0

import sys


directory = sys.argv[1]
file = "cluster_summary.tsv"
path = directory + "/" + file
ani = sys.argv[2]

dataset1_handler = open(path,"r")
fields_dict = {}
type_dict = {}
groups = []


lines_dt1 = dataset1_handler.readlines()


# Defining new cluster name
new_cluster_id = 1

print(directory, file=sys.stderr)

for l in lines_dt1[1:]:
    l = l.strip()
    fields = l.split("\t")

    del fields[2] # remove centroid column that I dont need

    if int(fields[2]) >= 75:
    	fields_dict[fields[1]] = fields
    	type_dict[fields[1]] = "NEW"
    	if fields[0] not in groups:
    		groups.append(fields[0])
    else:
    	fields_dict[fields[1]] = fields
    	type_dict[fields[1]] = "OUT"

if not groups:
	groups.append(1)

for keys in fields_dict.keys():
	if type_dict[keys] == "NEW":
 		print(*fields_dict[keys], 
 			directory + "-" + str(groups.index(fields_dict[keys][0])+1),
 			ani, "GOODBS", sep="\t")

out_group_number = str(groups.index(groups[-1])+2)

for keys in fields_dict.keys():
	if type_dict[keys] == "OUT":
		print(*fields_dict[keys], 
			directory + "-" + out_group_number + "-OUT",
			ani, "OUT", sep="\t")

# for i in range(len(groups)):
# 	print(i+1, groups[i])

dataset1_handler.close()
