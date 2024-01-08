DIR=$HOME/test
mkdir test

# get motifs directly from HOCOMOCOv11
cd $DIR
wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme 

# get hg38 FASTA by chromosome
FILE=$DIR/genome/chromFa.tar.gz
mkdir -p $(dirname $FILE)
wget -O $FILE http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar zxvf $FILE -C $(dirname $FILE)
echo $(seq 1 22) X Y M | tr ' ' '\n' | parallel -P1 cat $DIR/genome/chroms/chr{}.fa > $DIR/genome/hg38.fa

# Create scripts for each motif and chromosome
OUTDIR=$DIR/fimo_analysis
mkdir -p $OUTDIR/src/
OUT_SCRIPT=$OUTDIR/src/find_motifs.sh
rm -f ${OUT_SCRIPT}

FIMO=fimo
FASTA=$DIR/genome/hg38.fa

OUT=$OUTDIR/HOCOMOCOv11
mkdir -p $OUT
MOTIFS=$DIR/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
for MOTIF in $(grep MOTIF $MOTIFS | cut -f2 -d' ')
do
  mkdir -p $OUT/$MOTIF
  echo "$FIMO -o $OUT/$MOTIF/hg38 --motif $MOTIF --max-stored-scores 10000000 --thresh 1e-5 --no-qvalue $MOTIFS $FASTA" >> ${OUT_SCRIPT}
done

LOG=$OUTDIR/logs
mkdir -p $LOG $OUTDIR/jobs

cd $DIR
# convert GFF to starch for faster processing
rm -f $OUTDIR/src/convert_gff_to_starch.sh
rm -f $OUTDIR/src/convert_gff_to_bed.sh
for GFF in $(find $OUTDIR -name fimo.gff)
do
  STRCH=$(echo $GFF | sed 's/gff$/starch/g')
  BED=$(echo $GFF | sed 's/gff$/bed/g')
  echo "cat $GFF | gff2starch - > $STRCH" >> $OUTDIR/src/convert_gff_to_starch.sh
  echo "unstarch $STRCH > $BED" >> $OUTDIR/src/convert_gff_to_bed.sh
done

# convert to starch
cat $OUTDIR/src/convert_gff_to_starch.sh | parallel -P1

# convert startcfh to bed
cat $OUTDIR/src/convert_gff_to_bed.sh | parallel -P1
