# Installation

## Operating system

The program has been extensively tested and working on Linux systems. Due to the nature of basic dependencies, it's expected that it'll run on macOS and Windows systems as well with minor tweaks.

##  Software environment/ dependencies
 APP has been written in Perl. Perl comes pre-installed with most linux distributions. 
 
 Other dependencies include installation of the following programs. Individual setup scripts are provided for requirement-based installation.
  
 1. Command-line NCBI BLAST toolkit. If already installed, please ignore. To install: 
 ```
sh scripts/installation/setup_blast.sh
```
 2. Command-line NCBI Eutilities toolkit. If already installed, please ignore. To install: 
 ```
sh scripts/installation/setup_blast.sh
```
 
 3. Command-line TaxonKit NCBI Taxonomy toolkit. If already installed, please ignore. To install: 
 ```
sh scripts/installation/setup_blast.sh
```
 4. Command-line csvtk CSV/TSV toolkit. If already installed, please ignore. To install: 
 ```
sh scripts/installation/setup_blast.sh
```
 5. NCBI taxonomy database. If already installed, please ignore. To install: 
 ```
sh scripts/installation/setup_blast.sh
```


If none of the dependencies have been previously installed, please run: 

```
sh scripts/installation/complete_setup.sh
```
## Set permissions on file
```
sudo chmod 777 $HOME/bin_APP
```


