Installing MuDoGeR via conda can help the user to utilize only part of the workflow. However, it is more recommended for those with a deeper understanding of how conda environmets works, as manual adjustments may need to be made.

## Installation using Conda/Mamba
**1 - Install miniconda**

The MuDoGeR pipeline requires several tools that have multiple dependencies. More often than not, those tools require conflicting dependencies. To tackle this problem, MuDoGeR requires miniconda to be previously installed in your system. MuDoGeR will use miniconda to create multiple environments that are automatically orchestrated during the pipeline's functioning. Following you have a possible way of installing miniconda.

```console

#See documentation: https://docs.conda.io/en/latest/miniconda.html

$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ chmod +x Miniconda3-latest-Linux-x86_64.sh

$ ./Miniconda3-latest-Linux-x86_64.sh

$ export PATH=~/miniconda3/bin:$PATH

```

**2 - Install MuDoGeR**

Once you have miniconda installed and on your PATH, you can properly install MuDoGeR.
The MuDoGeR environment can be complex, and the installation script was designed to install and set up all necessary bioinformatics tools.

```console
#clone repository

$ git clone https://github.com/mdsufz/MuDoGeR.git

#Go to the MuDoGeR cloned repository folder
$ cd MuDoGeR

#Make sure you have conda ready and that you are in your base environment.
$ conda activate base
$ echo $CONDA_PREFIX

#You should see something like the following:
/path/to/miniconda3

#Install mamba in your base environment
$ conda install -c conda-forge mamba

#Run the installation script as follows
$ bash -i installation/install.sh

#Follow the instructions on the screen:
# Enter "y" if you want to install all modules, otherwise enter "n".
# If you entered "n",enter "y" for each of the modules you would like to install individually.

	The MuDoGeR's installation will begin..





	      (  )   (   )  )			
	       ) (   )  (  (			
	       ( )  (    ) )			
	       _____________			
	      <_____________> ___		
	      |             |/ _ \		
	      |               | | |		
	      |               |_| |		
	   ___|             |\___/		
	  /    \___________/    \		
	  \_____________________/		

	This might take a while. Time to grab a coffee...
```

**3 - Install necessary databases**

Several bioinformatics tools used within MuDoGeR require specific databases to work. We developed a database download and set up tool to make our lives easier. Make sure to run the database setup after MuDoGeR is installed.

You can also choose to install the databases used by each module individually. You can use the flag ```--dbs``` and choose the name of the module you want to set up the databases (all \[default], prokaryotes, viruses, eukaryotes).

Use this script if you want MuDoGeR to take care of everything. 

```console
#Make sure mudoger_env is activated. It should have been created when you ran 'bash -i installation/install.sh'
$ conda activate mudoger_env

#Go to MuDoGeR cloned directory
$ cd MuDoGeR

#Run database setup script
$ bash -i installation/database-setup.sh --dbs all -o /path/to/save/databases

#You can also check out the database-setup help information
$ bash -i installation/database-setup.sh --help

MuDoGeR database script v=1.0
Usage: bash -i database-setup.sh --dbs [module] -o output_folder_for_dbs

  --dbs all              		download and install the necessaries databases for all MuDoGeR modules [default]
  --dbs prokaryotes              	download and install the necessaries databases for prokaryotes module
  --dbs viruses              		download and install the necessaries databases for viruses module
  --dbs eukaryotes              	download and install the necessaries databases for eukaryotes module
  -o path/folder/to/save/dbs      	output folder where you want to save the downloaded databases
  --help | -h				show this help message
  --version | -v			show mudoger version


```

**3.1 - Update databases**

We plan to update the automatic database installation script at least once a year. The user can check the version of the database update script by using the ```--version``` flag. Finally, in case the user would like to manually update the databases, one should find the location of the databases used by MuDoGeR by running ```cat $CONDA_PREFIX/envs/mudoger_env/bin/database.sh```.

The user can then follow the instructions from the tool developer to update the desired database. The only requirement is that the root folder for the database uses the name of the tool as follows: ```buscodbs/  checkm/  checkv/  eukccdb/  gtdbtk/  vibrant/  wish/```.

**4 - Additional module 4 (eukaryotes) installation instructions**

Some tools used in module 4 (GENEMARK and MAKER2) require the user to provide information to the developers. Consequently, we could not implement an automatic installation and setup script. However, we created a tutorial to finish the module 4 setup.

The module 4 setup tutorial is found in [Module 4 setup](https://github.com/mdsufz/MuDoGeR/blob/master/installation/genemark_maker2_installation.md)
