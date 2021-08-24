# Installation

## Operating system

The program has been extensively tested and working on Linux systems. Due to the nature of basic dependencies, it's expected that it'll run on macOS and Windows systems as well with minor tweaks.

##  Software environment/ dependencies
 APP has been written in Perl. Perl comes pre-installed with most linux distributions. 
 
 Other dependencies include installation of the following programs. Individual setup scripts are provided for requirement-based installation.
  
 1. Command-line NCBI BLAST toolkit. If already installed, please ignore, else: 
 ```
sh scripts/installation/setup_blast.sh
```

 2. Command-line NCBI Eutilities toolkit. If already installed, please ignore, else: 
 ```
sh scripts/installation/setup_Eutilities.sh
```
 
 3. Command-line TaxonKit NCBI Taxonomy toolkit. If already installed, please ignore, else: 
 ```
sh scripts/installation/setup_TaxonKit.sh
```

 4. Command-line csvtk CSV/TSV toolkit. If already installed, please ignore, else: 
 ```
sh scripts/installation/setup_csvtk.sh
```

 5. NCBI taxonomy database. If already installed, please ignore, else: 
 ```
sh scripts/installation/setup_NCBI_taxDump.sh
```

If none of the dependencies have been previously installed, please run: 

```
sh scripts/installation/complete_setup.sh
```
## Perl dependencies

```
cpan App::cpanminus
cpanm List::MoreUtils 
```

## Set permission of files
```
sudo chmod 777 $HOME/bin_APP/*
```
## Set up taxonomic database for APP

The last step of installation involves construction of one-time taxonomic database compatible with APP. This database is to be created only the first time while running APP. With newer updates of NCBI database, this script can be run again to create an updated taxonomic file. 

Follow the steps as listed in [Taxonomy.md](https://github.com/sohamsg90/APP-Alieness-by-Phyletic-Pattern/blob/main/docs/Taxonomy.md) 
