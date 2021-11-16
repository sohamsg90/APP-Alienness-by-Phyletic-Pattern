##Create a directory for installing essentials

mkdir $HOME_PATH/bin_APP/

##Download the taxonkit executable 

wget -c -q https://github.com/shenwei356/taxonkit/releases/download/v0.8.0/taxonkit_linux_amd64.tar.gz
tar -zxvf taxonkit_linux_amd64.tar.gz
mv taxonkit /usr/local/bin
# cp taxonkit $HOME_PATH/bin_APP/

##Download the csvtk executable 

wget -c -q https://github.com/shenwei356/csvtk/releases/download/v0.23.0/csvtk_linux_amd64.tar.gz
tar -zxvf csvtk_linux_amd64.tar.gz
mv csvtk /usr/local/bin
#cp csvtk $HOME_PATH/bin_APP/

##Download the NCBI taxonomy database

wget -c -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz
mv names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME_PATH/bin_APP/


##Download the command-line BLAST toolkit

wget -c -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz
mv ncbi-blast-2.6.0+/bin/* /usr/local/bin
rm -R ncbi-blast-2.6.0+/

#wget -c ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh 

##Cpan modules
cpanm Data::Dumper
cpanm Getopt::Long
cpanm List::MoreUtils
cpanm Statistics::R

##Download and install edirect & dependencies
wget -c -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
#nquire -dwn ftp.ncbi.nlm.nih.gov entrez/entrezdirect edirect.tar.gz
tar -zxvf edirect.tar.gz
chmod +x edirect/*
mv edirect/* /usr/local/bin
rm -R edirect

wget -c -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz
#nquire -dwn ftp.ncbi.nlm.nih.gov entrez/entrezdirect xtract.Linux.gz
gunzip -f xtract.Linux.gz
chmod +x xtract.Linux
mv xtract.Linux xtract 
mv xtract /usr/local/bin/

wget -c -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/rchive.Linux.gz
#nquire -dwn ftp.ncbi.nlm.nih.gov entrez/entrezdirect rchive.Linux.gz
gunzip -f rchive.Linux.gz
chmod +x rchive.Linux
mv rchive.Linux rchive 
mv rchive /usr/local/bin/

wget -c -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/transmute.Linux.gz
#nquire -dwn ftp.ncbi.nlm.nih.gov entrez/entrezdirect transmute.Linux.gz
gunzip -f transmute.Linux.gz
chmod +x transmute.Linux
mv transmute.Linux transmute 
mv transmute /usr/local/bin/

##Download and install cgviewer
wget -c -q https://github.com/paulstothard/cgview/archive/refs/tags/v2.0.2.tar.gz
tar -zxvf v2.0.2.tar.gz
chmod +x cgview-2.0.2/*
cp cgview-2.0.2/bin/cgview.jar /usr/local/bin/
cp cgview-2.0.2/scripts/cgview_xml_builder/cgview_xml_builder.pl /usr/local/bin/
rm -R cgview-2.0.2

##Copy APP programs to /usr/local/bin
cd APP_scripts/
chmod +x *
cp * /usr/local/bin
cd ..
##Remove compressed files
rm *.gz *.dmp *.sh gc.prt *.txt


