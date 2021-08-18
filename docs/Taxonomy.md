# Setting up the taxonomy database
A pre-built taxonomy database is provided along with the main program. The database has been updated as of Aug-15-2021. 

## Steps to update the pre-built taxonomy database
A bash script has been included to download the Refseq assembly details from NCBI FTP-server. Post-processing using linux commands, and, by means of TaxonKit, the database can be compiled. The following command can be used to generate the same

`sh setup_taxonomy_database.sh`

The output file needs to be copied and pasted to the same directory housing the main script `APP.pl`. Discard the old database file prior to updating the database.


