#!/bin/sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=solomonb@stanford.edu
#SBATCH --time=13-23:05 # Runtime in D-HH:MM
#SBATCH --job-name=SCHLA
#SBATCH --nodes=1 # Ensure that all cores are reserved on one machine
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7G
#SBATCH --partition=khatrilab # Partition allocated for the lab
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# SET GLOBAL VARIABLES
# General
export GIT_DIR=/labs/khatrilab/solomonb/hla_project/hla_benchmark/data/isb/scHLAcount
export BASE_DIR=/labs/khatrilab/solomonb/covid/isb/scHLAcount
export LOG_DIR=$GIT_DIR/logs/$(date +'%y%m%d_%H%M%S')
# FASTQ/HISAT
export INDEX_DIR=/labs/khatrilab/solomonb/rnaseq_processing/hisat2/hisat_arcas/hisat_data/grch38
export BAM_DIR=$BASE_DIR/bam
# HLA references
export HLA_DIR=$GIT_DIR/hla_references
export HLANUC=/labs/khatrilab/solomonb/references/IMGTHLA/hla_nuc.fasta
export HLAGEN=/labs/khatrilab/solomonb/references/IMGTHLA/hla_gen.fasta
# SCHLA
export BARCODE_DIR=$GIT_DIR/barcodes
export GENOTYPE_DIR=$GIT_DIR/genotypes
export SCHLACOUNT_DIR=$GIT_DIR/output
export TEMP_DIR=$GIT_DIR/temp_fastq
# SLURM 
export N_CORES=$SLURM_CPUS_PER_TASK


# CREATE DIRECTORIES
if [ ! -d $LOG_DIR ]; then mkdir -p $LOG_DIR;fi
if [ ! -d $BAM_DIR ]; then mkdir -p $BAM_DIR;fi
if [ ! -d $SCHLACOUNT_DIR ]; then mkdir -p $SCHLACOUNT_DIR;fi
if [ ! -d $TEMP_DIR ]; then mkdir -p $TEMP_DIR;fi

# CREATE HLA REFERECE ##########################################################
HLA_REFERENCE(){
  source /labs/khatrilab/solomonb/miniconda3/etc/profile.d/conda.sh
  conda activate samtools
  
  printf "\n\
########################## CREATE REFERENCE ####################################\
  \n" >> $LOG_DIR/${1}/${1}_${2}.log
  printf "\n### [START___REFERENCE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_${2}.log
  
  [ ! -d $HLA_DIR/${2} ] && mkdir -p $HLA_DIR/${2}
  
  while read -r line; do grep -F -m 1 $line $HLANUC >> $HLA_DIR/${2}/${1}_tmpallele.txt; done < $GENOTYPE_DIR/${2}/${1}_hla.tsv
  
  samtools faidx  $HLANUC $(cut -f1 -d' ' $HLA_DIR/${2}/${1}_tmpallele.txt | tr '>' ' ' | tr '\n' ' ') > $HLA_DIR/${2}/${1}_cds.fasta 2>> $LOG_DIR/${1}/${1}_${2}.log
  samtools faidx  $HLAGEN $(cut -f1 -d' ' $HLA_DIR/${2}/${1}_tmpallele.txt | tr '>' ' ' | tr '\n' ' ') > $HLA_DIR/${2}/${1}_gen.fasta 2>> $LOG_DIR/${1}/${1}_${2}.log
  
  while read -r line; do IFS=' '; read -r f1 f2 <<<"$line"; sed -i"" "s/$f1/$f1 $f2/g" $HLA_DIR/${2}/${1}_cds.fasta; done < $HLA_DIR/${2}/${1}_tmpallele.txt 2>> $LOG_DIR/${1}/${1}_${2}.log
  while read -r line; do IFS=' '; read -r f1 f2 f3 <<<"$line"; sed -i"" "s/$f1/$f1 $f2/g" $HLA_DIR/${2}/${1}_gen.fasta; done < $HLA_DIR/${2}/${1}_tmpallele.txt 2>> $LOG_DIR/${1}/${1}_${2}.log
  
  rm $HLA_DIR/${2}/${1}_tmpallele.txt
  printf "### [COMPLETE___REFERENCE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_${2}.log
}
export -f HLA_REFERENCE


# DEFINE scHLA GENOTYPING PIPELINE #############################################
SCHLACOUNT(){
  source /labs/khatrilab/solomonb/miniconda3/etc/profile.d/conda.sh
  conda activate samtools
  
  printf "\n\
########################## RUN SCHLACOUNT ####################################\
  \n" >> $LOG_DIR/${1}/${1}_${2}.log
  printf "\n### [START___SCHLACOUNT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_${2}.log
  
  echo "### Starting  scHLAcount at $(date +'%D %X')" >> $LOG_DIR/${1}/${1}_${2}.log
  
  # if [ -d $SCHLACOUNT_DIR/${1}_results ]; then rm -r $SCHLACOUNT_DIR/${1}_results;fi
  [ ! -d SCHLACOUNT_DIR/${2} ] && mkdir -p $SCHLACOUNT_DIR/${2}

  sc_hla_count \
  --bam $BAM_DIR/${1}.bam \
  --cell-barcodes $BARCODE_DIR/${1}_barcode.tsv \
  --out-dir $SCHLACOUNT_DIR/${2}/${1}_results \
  --fasta-cds $HLA_DIR/${2}/${1}_cds.fasta \
  --fasta-genomic $HLA_DIR/${2}/${1}_gen.fasta\
  >> $LOG_DIR/${1}/${1}_${2}.log \
  2>> $LOG_DIR/${1}/${1}_${2}.log

  printf "### [COMPLETE___SCHLACOUNT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_${2}.log
}
export -f SCHLACOUNT


# DEFINE PIPELINE CONTROLLER FUNCTION ##########################################
PIPELINE(){
  echo "START: sample $1 at $(date +'%D %X')"
  for G in arcasHLA hlaminer invitro optitype phlat
  do
    [ ! -d $LOG_DIR/${1} ] && mkdir -p $LOG_DIR/${1}
    HLA_REFERENCE $1 $G
    SCHLACOUNT $1 $G
  done
  echo "COMPLETE: sample $1 at $(date +'%D %X')"
}
export -f PIPELINE

cat $GIT_DIR/BL_fastq_files.txt | parallel --delay 15 -j $SLURM_NTASKS --joblog $LOG_DIR/parallel.log PIPELINE {}