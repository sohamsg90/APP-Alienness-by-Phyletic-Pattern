##Use the sampe multiFasta file in example folder to run a sample run

#Example taxonomic file provided here. If not, use: 
#perl scripts/taxonomy/input_file_prepare.pl example/accession_number.txt

#Then run
cp database/1.REFSEQ_all_kingdoms_ftplinks_lineage.txt .
perl APP.pl -q example/NC_004088.fa -t example/3.query_speciesID_taxID.txt -o NC_004088 -f multifasta

