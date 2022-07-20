# Finish Module 4 installation

### Activate MuDoGeR environment and load config file

```console
$ conda activate mudoger_env
$ config_file="$(which config.sh)"
$ source "$config_file"

```
## Complete GENEMARK installation

1. ACCESS GENEMARK WEBPAGE

  http://opal.biology.gatech.edu/GeneMark/license_download.cgi


2. SELECT OPTIONS 

      **GeneMark-ES/ET/EP ver \*_lic and LINUX 64**


3. FILL IN THE CREDENTIALS WITH ***YOUR NAME, E-MAIL, INSTITUTION, ETC... ***


4. CLICK ON ***'I agree the terms of this license agreement'***


5. DOWNLOAD THE 64_bit key and the program files provided

      It should look something like the following:

        ```console
        $ wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_HZzc0/gmes_linux_64.tar.gz
        $ wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_HZzc0/gm_key_64.gz
        ```

      You should have the following files:
      **gmes_linux_64.tar.gz** and  **gm_key_64.gz**


6. RETRIEVE PATH TO GENEMARK ENV (this environment is already created and ready once you installed module 4 tools by running install.sh. See [Installation](https://github.com/JotaKas/MuDoGeR/blob/master/README.md#installation))

```console
$ genemarker_environment="$MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env"
```

7. COPY THE FILES TO THE FOLLOWING SPECIFIC FOLDER (remember to first run the install.sh script)

```console
$ cp /path/to/gmes_linux_64.tar.gz  $genemarker_environment/GENEMARK_MAIN
$ cp /path/to/gm_key_64.gz  $genemarker_environment/GENEMARK_MAIN

# IT SHOULD LOOK LIKE THIS:
$ ls  $genemarker_environment/GENEMARK_MAIN
$ gmes_linux_64.tar.gz  gm_key_64.gz
```

8. ENTER THE FOLDER AND DECOMPRESS FILES

```console
$ cd  $genemarker_environment/GENEMARK_MAIN
$ gunzip *
```

9. COPY AND RENAME KEY FILE TO YOUR HOME FOLDER 

```console
$ cp gm_key_64 ~/.gm_key
```

10. DECOMPRESS TOOL

```console
$ tar -xvf gmes_linux_64.tar
$ rm -fr gmes_linux_64.tar # Make sure you remove the tar file once extracted
```

11. ACTIVATE THE genemark_env CONDA ENVIRONMENT (created during the installation procedure)

```console
$ conda activate $genemarker_environment
```

12. GET PERL PATH

```console
$ perl_path="$(which perl)"
$ cd - #return
```

13. CONFIGURE GENEMARK TO USE THE CORRECT PERL

```console
$ cd $genemarker_environment/GENEMARK_MAIN/gmes_linux_64
$ perl change_path_in_perl_scripts.pl "$perl_path"
$ cd - #return
```

14. Return to mudoger_env
```console
$ conda deactivate
```


## COMPLETING INSTALLATION FOR MAKER2

1. DOWNLOAD MAKER (Registration required)

    Go to http://weatherby.genetics.utah.edu/cgi-bin/registration/maker_license.cgi and fill out your information to download the tool
    Copy the download link provided (Usually by Right-click -> 'Copy Link') and type:
    ```console
    $ wget *Download link*

    #Should be something like:
    # $ wget http://weatherby.genetics.utah.edu/maker_downloads/1BF1/66D1/3FAA/FF20FB345D6221AA8A338D1B9D8A/maker-3.01.04.tgz
    ```


2. RETRIEVE PATH TO MUDOGER_ENV (this environment is already created and ready once you used install.sh. See [Installation](https://github.com/JotaKas/MuDoGeR/blob/master/README.md#installation))

```console
$ maker2_environment="$MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env" 
```


3. PLACE THE DOWNLOADED MAKER FILE INSIDE SPECIFIC FOLDER created during install.sh (remember to first run the install.sh script)
```console
$ cp maker-3.01.04.tgz $maker2_environment/MAKER_MAIN
```


4. ENTER FOLDER AND DECOMPRESS FILES
```console
$ cd  $maker2_environment/MAKER_MAIN
$ tar xvfz maker-3.01.04.tgz
```


5. CREATE A EXE FOLDER INSIDE THE MAKER EXTRACTED FOLDER
```console
$ mkdir $maker2_environment/MAKER_MAIN/maker/exe
```


6. DOWNLOAD DFAM DATABASE. You should have enough space for this database where conda is installed (~97GB)
```console
$ wget -P $maker2_environment/MAKER_MAIN/ https://www.dfam.org/releases/Dfam_3.5/families/Dfam.h5.gz
$ gunzip Dfam.h5.gz
```


7. DOWNLOAD AND MOVE SNAP
```console
$ wget -P $maker2_environment/MAKER_MAIN/ http://korflab.ucdavis.edu/Software/snap-2013-11-29.tar.gz
$ tar xvfz snap-2013-11-29.tar.gz
$ mv snap $maker2_environment/MAKER_MAIN/maker/exe
```


8. COPY SOME SCRIPTS FROM YOUR GENEMARK INSTALLATION (Installation done [previously](https://github.com/JotaKas/MuDoGeR/blob/master/installation/genemark_maker2_installation.md#complete-genemark-installation), just copied the necessaries bin)
```console
$ cp -r $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env/GENEMARK_MAIN/gmes_linux_64/* $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env/bin/
```


9. DOWNLOAD REPEATMASKER
```console
$ wget -P $maker2_environment/MAKER_MAIN/ https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz
```


10. EXTRACT REPEATMASKER
```console
$ tar xvfz RepeatMasker-4.1.2-p1.tar.gz
```


11. MOVE DFAM DATABASE TO LIBRARIES FOLDER IN REPEATMASKER
```console
$ mv -f $maker2_environment/MAKER_MAIN/Dfam.h5 $maker2_environment/MAKER_MAIN/RepeatMasker/Libraries/
```

12. MOVE REPEATMASKER TO EXE FOLDER IN MAKER
```console
$ mv $maker2_environment/MAKER_MAIN/RepeatMasker $maker2_environment/MAKER_MAIN/maker/exe/
$ cd - #return
```


13. INSTALL REPEATMASKER
```console
$ conda activate $maker2_environment
$ cd $maker2_environment/MAKER_MAIN/maker/exe/RepeatMasker/
```

   You will be prompt to enter some tools paths to configure repeatmasker.
   Maker sure your maker2 enviroment is activated. It was created during initial [installation](https://github.com/JotaKas/MuDoGeR/blob/master/README.md#installation)
    Ideally it should automatically find everything. If that is the case, just hit enter when prompt.
    
```console
$ perl ./configure
```
   Hit enter if TRF path is correct
   Select 3 (HMMER) as search engine and set as default (Y)
   Hit 5 to be done and wait database configuration
```console
$ cd - #return
```


14. CONFIGURE MAKER
```console
$ cd $maker2_environment/MAKER_MAIN/maker/src
$ perl Build.PL #(Currently answer N to set MPI)
```
   If any Perl dependency was not properly installed, for any reason, please run the following for the each dependency
```console
$ cpanm FAILED::DEPENDENCY
$ #ex: cpanm Bit::Vector
```
  If that did not work, for any reason, use the maker2 solution and configure maker2 again
  ```console
  ./Build installdeps
  ```


15. INSTALL MAKER
```console
$ ./Build install

$ export ZOE="$maker2_environment/MAKER_MAIN/maker/exe/snap/Zoe"
$ export AUGUSTUS_CONFIG_PATH="$maker2_environment/config"
```


16. ADD MAKER TO PATH
```console
$ export PATH=$maker2_environment/MAKER_MAIN/maker/bin:$PATH
$ cd - #return
```


17. Return to mudoger_env
```console
$ conda deactivate
```
