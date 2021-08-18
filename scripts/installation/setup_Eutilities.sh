##Download the NCBI Eutilities toolkit. 
##Detailed manual can be accessed From: https://www.ncbi.nlm.nih.gov/books/NBK179288/
cd ~
/bin/bash
perl -MNet::FTP -e \
'$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
$ftp->login; $ftp->binary;
$ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
gunzip -c edirect.tar.gz | tar xf -
rm edirect.tar.gz
builtin exit
mkdir -p $HOME/bin_APP/edirect
export PATH=${PATH}:$HOME/bin_APP/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/bin_APP/edirect"
./edirect/setup.sh
