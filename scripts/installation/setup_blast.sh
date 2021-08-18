##Download the command-line BLAST toolkit

wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz
mkdir -p $HOME/bin_APP/blast
cp ncbi-blast-2.6.0+/bin/* $HOME/bin_APP/blast
