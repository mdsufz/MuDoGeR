

mudoger_lettering(){
echo -e "\n"
echo -e "\t███    ███ ██    ██ ██████   ██████   ██████  ███████ ██████  "
echo -e "\t████  ████ ██    ██ ██   ██ ██    ██ ██       ██      ██   ██ "
echo -e "\t██ ████ ██ ██    ██ ██   ██ ██    ██ ██   ███ █████   ██████  "
echo -e "\t██  ██  ██ ██    ██ ██   ██ ██    ██ ██    ██ ██      ██   ██ "
echo -e "\t██      ██  ██████  ██████   ██████   ██████  ███████ ██   ██ "                                                                                                                    
echo -e "\t\t\tMulti-Domain Genome Recovery"
echo -e "\t\t\t\tVersion 1.0.0\n\n"
}
  

mudoger_lettering

echo '--> Please provide a path to the metadata file. ($ mudoger -help -md) for explanation'

read metadata
echo 'path to metadata file: '$metadata
echo 'reading '$metadata

out="$(python MuDoGeR/tools/mdcheck.py $metadata)"

echo -e "$out"

### check if there was any problem with metadata file
if [[ $out =~ "Closing" ]]; then
   echo -e "\n TIP: run \"mudoger -h\" for help "
   exit 0
fi


# if everything is okay, proceeds with asking user for output folder
echo '--> Please provide a path for the output folder. If it does not exist, it will be created'
read output_folder

mkdir -p "$output_folder"
echo '--> Folder created'
