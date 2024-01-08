#!/bin/bash

mem="16GB"
time="12:00:00"
partition="shortq"
ncpus=8

DIR=resources/TFBS_5e-4
mkdir -p $DIR
mkdir -p $DIR/src/
LOG_DIR=$DIR/logs
mkdir -p ${LOG_DIR}

# wget -O $GENOME http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
# tar zxvf $GENOME -C $(dirname $GENOME)
# echo $(seq 1 22) X Y M | tr ' ' '\n' | parallel -P1 cat $DIR/genome/chroms/chr{}.fa > $GENOME
GENOME=$DIR/genomes/hg38.fa

for collection in "HOCOMOCOv11_core_HUMAN_mono" # "JASPAR2022_CORE_vertebrates_non-redundant"
do
  MOTIFS=$DIR/$collection.meme
  OUT=$DIR/$collection/hg38
  mkdir -p $OUT $OUT/regions
  OUT_SCRIPT=$DIR/src/${collection}.sh
  rm -f ${OUT_SCRIPT}
  
  for motif in $(grep MOTIF $MOTIFS | cut -f2 -d' ') 
  do
    echo "fimo -oc $OUT/$motif --motif $motif --max-stored-scores 50000000 --thresh 5e-4 --no-qvalue $MOTIFS $GENOME; cat $OUT/$motif/fimo.gff | gff2starch - > $OUT/$motif/$motif.starch; unstarch $OUT/$motif/$motif.starch | cut -f 1-3 > $OUT/regions/$motif.bed; bedtools merge -i $OUT/regions/$motif.bed > $OUT/regions/$motif.tmp; mv $OUT/regions/$motif.tmp $OUT/regions/$motif.bed" >> ${OUT_SCRIPT}
  done
  
  cat ${OUT_SCRIPT} | while read cmd; 
  do  
	  job=$(echo $cmd | cut -f5 -d' ')
	  echo "Submitting $job"
	  echo $cmd
	  sbatch --wrap "${cmd}" --job-name="FIMO_${job}" --output="${LOG_DIR}/${job}.log" \
    --mem=${mem} --cpus-per-task=${ncpus} --time=${time} --partition=${partition} --qos=${partition}
  done
done
