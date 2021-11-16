#!/bin/bash
FILE=$1

JoinIntoGroupsOf() {
  xargs -n "$@" echo |
  sed 's/ /,/g'
}

cat $FILE | JoinIntoGroupsOf 5000 | xargs -n 1 sh -c \
 'esummary -db nuccore -id $0 | xtract -pattern DocumentSummary -element Caption,TaxId ' \
 > $FILE\_taxID.txt

awk -F "\t" '{print $2}' $FILE\_taxID.txt > tax_ID_of_interest.txt
rm $FILE\_taxID.txt

taxonkit lineage tax_ID_of_interest.txt --data-dir $HOME_PATH/bin_APP/ > tax_ID_of_interest_lineage.txt

##For full names
cat tax_ID_of_interest_lineage.txt \
    | taxonkit reformat --data-dir $HOME_PATH/bin_APP/ -F  \
    | csvtk -H -t cut -f 1,3 \
    | csvtk -H -t sep -f 2 -s ';' -R \
    | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species \
    | csvtk pretty -T > tax_ID_of_interest_lineage_ID_format.txt 

##For numeric ID format
cat tax_ID_of_interest_lineage.txt  \
    | taxonkit reformat --data-dir $HOME_PATH/bin_APP/ -t  \
    | csvtk -H -t cut -f 1,4 \
    | csvtk -H -t sep -f 2 -s ';' -R \
    | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species \
    | csvtk pretty -T > tax_ID_of_interest_lineage_ID_format.txt 


