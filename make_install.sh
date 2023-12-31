mkdir tool
mkdir reference

# Download blast from ncbi
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz -P tool/

# Install blastn
cd tool
tar -xf ncbi-blast-2.6.0+-x64-linux.tar.gz
mv ncbi-blast-2.6.0+/bin/blastn .
rm -rf ncbi-blast-2.6.0+*
cd ..

# Download genome fasta.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz -P reference/
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -P reference/

# Decompress fa.gz
gunzip reference/hg19.fa.gz
gunzip reference/hg38.fa.gz

#rm reference/*.gz
