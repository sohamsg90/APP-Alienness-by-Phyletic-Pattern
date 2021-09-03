##Create a directory for installing essentials

mkdir $HOME/bin_APP/

##Download the taxonkit executable 

wget -c -v https://github.com/shenwei356/taxonkit/releases/download/v0.8.0/taxonkit_linux_amd64.tar.gz
tar -zxvf taxonkit_linux_amd64.tar.gz
cp taxonkit $HOME/bin_APP/

##Download the csvtk executable 

wget -c -v https://github.com/shenwei356/csvtk/releases/download/v0.23.0/csvtk_linux_amd64.tar.gz
tar -zxvf csvtk_linux_amd64.tar.gz
cp csvtk $HOME/bin_APP/

##Download the NCBI taxonomy database

wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/bin_APP/





##Download the command-line BLAST toolkit

wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz
mkdir -p $HOME/bin_APP/blast
cp ncbi-blast-2.6.0+/bin/* $HOME/bin_APP/blast

# ##Download the NCBI Eutilities toolkit. From: https://www.ncbi.nlm.nih.gov/books/NBK179288/
echo "For installation of NCBI Eutilities, please read and execute the code written in the file setup_Eutilites.txt "
chmod 777 $HOME/bin_APP/*
export PATH="$HOME/bin_APP:$PATH"
