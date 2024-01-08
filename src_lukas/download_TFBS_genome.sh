DIR=resources/TFBS
mkdir $DIR

# get hg38 FASTA by chromosome
FILE=$DIR/genome/chromFa.tar.gz
mkdir -p $(dirname $FILE)
wget -O $FILE http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar zxvf $FILE -C $(dirname $FILE)

# combine chromosomes
echo $(seq 1 22) X Y M | tr ' ' '\n' | parallel -P1 cat $DIR/genome/chroms/chr{}.fa > $DIR/genome/hg38.fa

