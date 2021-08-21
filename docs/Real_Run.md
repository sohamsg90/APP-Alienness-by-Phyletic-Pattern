# Describe how to run program

## Perl dependencies
 List::MoreUtils (https://metacpan.org/pod/List::MoreUtils)

## Sample run
The associated dependencies of APP must be installed prior to attempt running the program. The details of these have been provided in [Installation.md](https://github.com/sohamsg90/APP-Alieness-by-Phyletic-Pattern/blob/main/docs/Installation.md). 

## Input
APP requires - 
1. A multiFASTA protein sequence file and 
2. A file with taxonomic IDs of the associated genome, as the primary input files. 

The taxonomic IDs can be obtained by using the script `input_file_prepare.pl`. An input file with the genome accession number (NCBI) is to be provided like,

`perl scripts/taxonomy/input_file_prepare.pl accession_number.txt`

This generates a file with the species ID and taxonomic ID of the genome (from NCBI taxonomy database).

The output of the file looks like this:







