#### LITTLE HELP ####
## All conda environments is created inside
## dependencies/conda/envs directory.
## You can manually activate the environment
## using the following command
## $ conda activate path/to/mudoger_env/dependencies/conda/envs/environment_name

#### LITTLE HELP ####
## All modules installed here will be named as "name_env"
## Ex.: metawrap_env

################# INSTALLATION'S PRE-CONFIGURATION #################

## Here is where MuDoGeR's pre-configuration will 
## be defined. In general, in config.sh we collect 
## some MuDoGeR's important paths

## Importing config.sh
#echo dirname "$dirname"
#echo '----'

#conda activate base
mamba install -y -c anaconda conda-package-handling
mamba install -y -c conda-forge gawk 
mamba install -y libarchive==3.5.2 -c conda-forge

# Check if conda is installed
if command -v conda > /dev/null; then
  echo "Conda is installed."
else
  echo "Conda is not installed. Please install Conda and try again."
  exit 1
fi

##Required to be installed:
# mamba in base environment

# Check if mamba is installed
if conda list | grep -q mamba; then
  echo "Mamba is already installed. Version:"
  mamba --version
else
  echo "Mamba is not installed. Please install it before proceeding"
  echo "For instance, you can use:"
  echo "conda install -c conda-forge mamba"
  echo "while you have your base environment activated"
  exit 1
fi


MUDOGER_CONDA_ENVIRONMENT_PATH=$CONDA_PREFIX/envs/mudoger_env
echo MUDOGER_CONDA_ENVIRONMENT_PATH="$MUDOGER_CONDA_ENVIRONMENT_PATH" > $(dirname $0)/temp
cat $(dirname $0)/temp  $(dirname $0)/.config_std.sh > $(dirname $0)/config.sh
rm -f $(dirname $0)/temp 


source $(dirname $0)/config.sh
source $(dirname $0)/installation_utils.sh


## Moving the installation scripts to MuDoGeR's environment
# check if mudoger_env already exists
verify_if_mudoger_env_exists "$MUDOGER_CONDA_ENVIRONMENT_PATH"
if [ $mudoger_there == "yes" ]  # if yes, skip installation
then  
echo "Skip mudoger conda env installation. Moving forward..."
else 
echo "Installing mudoger conda env..."  # if no, move forward
start_pre_configuration
fi


echo "your MuDoGeR's path is $MUDOGER_CONDA_ENVIRONMENT_PATH"


yes | cp -rf $DEPENDENCIES_SCRIPTS_PATH $MUDOGER_CONDA_ENVIRONMENT_PATH

yes | cp -rf $INSTALLATION_SCRIPTS_PATH $MUDOGER_DEPENDENCIES_PATH


conda activate mudoger_env
core_env="$(echo $PATH | cut -f1 -d':')"


source installation/config.sh
source installation/installation_utils.sh


chmod +x $MUDOGER_WORK_SCRIPTS/*
chmod +x $MUDOGER_WORK_MODULES/*
chmod +x $MUDOGER_MASTER_SCRIPT 
chmod +x $(dirname $0)/config.sh 

yes | cp -rf $MUDOGER_WORK_SCRIPTS/*  $MUDOGER_CONDA_ENVIRONMENT_PATH/bin
yes | cp -rf $MUDOGER_WORK_MODULES/*  $MUDOGER_CONDA_ENVIRONMENT_PATH/bin
yes | cp -rf $MUDOGER_MASTER_SCRIPT   $MUDOGER_CONDA_ENVIRONMENT_PATH/bin
yes | cp -rf $(dirname $0)/config.sh  $MUDOGER_CONDA_ENVIRONMENT_PATH/bin

################# CHOOSING WICH MODULE TO INSTALL #################

## Giving the user the option of which modules he wants to install
## The Module 1. Pre-processing is the base of the other modules in MuDoGeR
## it's installation is automatically set

echo -e "\n### WELCOME TO MuDoGeR! ###\n"
echo "Do you want to install all MoDuGeR's Acessories Modules?"
echo "- Module 1. Pre-Processing"
echo "- Module 2. Recovery of Prokaryotic MAGs"
echo "- Module 3. Uncultivated viral genomes"
echo "- Module 4. Eukaryotic MAGs"
echo "- Module 5. Relative abundance and genome coverage Table"
while :
do
	read choose
	if [ $choose = y -o $choose = Y ];
	then
            	echo "------> Installing all modules"
		install_module_1_option=$choose
		install_module_2_option=$choose
		install_module_3_option=$choose
		install_module_4_option=$choose
		install_module_5_option=$choose
		break
	elif [ $choose = n -o $choose = N ]
	then
		echo "Choose the Modules you want to install"
		break
	else
		echo "Command not found, please, try again"
	fi
done
if [ $choose = n -o $choose = N ];
then
	echo "Do you want to install Module 1 (Pre-Processing)? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_module_1_option=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Module 1 denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done
	
	echo "Do you want to install Module 2 (Recovery of Prokaryotic MAGs)? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_module_2_option=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Module 2 denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done

	echo "Do you want to install Module 3 (Uncultivated viral genomes)? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_module_3_option=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Module 3 denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done

	echo "Do you want to install Module 4 (Eukaryotic MAGs)? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_module_4_option=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Module 4 denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done

	echo "Do you want to install Module 5 (Relative abundance and genome coverage Table)? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_module_5_option=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Module 5 denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done
fi

################# INSTALLING CHOSEN MODULES #################

## The Module 1. Pre-processing is the base of the other modules in MuDoGeR
## it's installation is automatically set

echo -e "\nThe MuDoGeR's installation will begin..\n"


coffe_time


if [ ! -z $install_module_1_option ];
then
	echo "-----> installing module 1"
	call_installation_script install_module_1
fi

if [ ! -z $install_module_2_option ];
then
	echo "-----> installing module 2"
	call_installation_script install_module_2
fi
if [ ! -z $install_module_3_option ];
then
    	echo "-----> installing module 3"
	call_installation_script install_module_3
fi
if [ ! -z $install_module_4_option ];
then
    	echo "-----> installing module 4"
	call_installation_script install_module_4
fi
if [ ! -z $install_module_5_option ];
then
    	echo "-----> installing module 5"
	call_installation_script install_module_5
fi
