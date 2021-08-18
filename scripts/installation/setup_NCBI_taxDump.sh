wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/bin_APP/
