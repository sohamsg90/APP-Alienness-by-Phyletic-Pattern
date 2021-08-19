##Download assembly report of refseq from NCBI ftp-server
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

##Extract recods which are either "Complete Genome" or "Chromosome"
awk -F "\t" '$12=="Complete Genome" || $12=="Chromosome" && $11=="latest"{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$20}' assembly_summary_refseq.txt > ftpdirpaths.txt

##Extract taxids from records in order. Not sorted or duplicates removed, to maintain order.
awk -F "\t" '{print $1}' ftpdirpaths.txt >  all_RefSeq_sequential_TAXIDS.txt

##Extract species name and ftp links for later
awk -F "\t" '{print $1"\t"$3"\t"$8}' ftpdirpaths.txt > extracted_ftpdirpaths.txt
#Add headers to file
sed -i '1i taxid\tspecies_name\tftp_links' extracted_ftpdirpaths.txt

##Extract lineage using TaxonKit 
taxonkit lineage all_RefSeq_sequential_TAXIDS.txt --data-dir $HOME/bin_APP/ > all_RefSeq_sequential_TAXIDS_lineage.txt

##Extract taxonomic id lineage using TaxonKit
cat all_RefSeq_sequential_TAXIDS_lineage.txt  \
    | taxonkit reformat --data-dir $HOME/bin_APP/ -t  \
    | csvtk -H -t cut -f 1,4 \
    | csvtk -H -t sep -f 2 -s ';' -R \
    | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species \
    | csvtk pretty -T > all_RefSeq_sequential_TAXIDS_lineage_byID.txt
    
##Extract full name lineage using TaxonKit
cat all_RefSeq_sequential_TAXIDS_lineage.txt  \
    | taxonkit reformat --data-dir $HOME/bin_APP/ -t  \
    | csvtk -H -t cut -f 1,3 \
    | csvtk -H -t sep -f 2 -s ';' -R \
    | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species \
    | csvtk pretty -T > all_RefSeq_sequential_TAXIDS_lineage_byWord.txt
    
##Combine all files together
paste -d "\t" all_RefSeq_sequential_TAXIDS_lineage_byWord.txt all_RefSeq_sequential_TAXIDS_lineage_byID.txt extracted_ftpdirpaths.txt > temp1.txt
##Remove repeatative columns and generate final file
awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19}' temp1.txt >  1.REFSEQ_all_kingdoms_ftplinks_lineage.txt

##Clean up
rm all_* assembly* extracted* ftpdir* taxIDs* temp1.txt

