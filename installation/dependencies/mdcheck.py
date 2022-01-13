
import sys
import os

md_file=sys.argv[1]
#md_file='/home/masun/Documents/MDS/projects/mdt/scripts/metadata_example-1.csv'
ids_dic,files,groups_dic={},{},{}
f=open(md_file,'r')

if '.csv' in md_file:
    sep=','
if '.tsv'  in md_file:
    sep='\t'

exit=0
while True:
    l=f.readline()
    if not l: break
    else:
       # print l
        id,filepath,group=l.strip().split(sep)

        isfile=os.path.isfile(filepath.replace("\"",""))
        if isfile == False:
            print '(\n--> Oops, it looks like the metadata file provided is incorrect.)'
            print '(    Please make sure all file paths are correct. Closing now.)'
            exit=1
            break

        else:
            if id not in ids_dic:
                ids_dic[id]=[filepath]
            else:
                ids_dic[id].append(filepath)

            if group not in groups_dic:
                groups_dic[group]=[filepath]
            else:
                groups_dic[group].append(filepath)

f.close()

if exit != 1:
    print '(\n-> Files are present')
    print '(--> Number of samples',len(ids_dic))
    print '(--> Number of groups',len(groups_dic))
