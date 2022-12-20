#!/bin/sh
#SBATCH --job-name=HLA_GEN
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7G

# SET GLOBAL VARIABLES
### Base directory that contains a directory with raw sequence files
### Additional data and log directories will be created in this directory
export BASE_DIR=
### Name of the directory within base_dir containing raw sequence files. Currently `raw_fastq`
export SAMPLE_DIR=$BASE_DIR/raw_fastq
### Name of file containing a list of samples IDs to be processed
export SAMPLE_LIST=$BASE_DIR/all_samples.txt
### Names for output directories that will be greated
export LOG_DIR=$BASE_DIR/logs/$(date +'%y%m%d_%H%M%S')
export TEMP_DIR=$BASE_DIR/temp_fastq
export FASTQC_DIR=$BASE_DIR/fastqc
export BAM_DIR=$BASE_DIR/bam
export ARCAS_DIR=$BASE_DIR/arcasHLA
export PHLAT_DIR=$BASE_DIR/phlat
export OPTITYPE_DIR=$BASE_DIR/optitype
export HLAMINER_DIR=$BASE_DIR/hla_miner

### Location of software directories for Bowtie, PHLAT, OptiType, etc..
export SOFTWARE_DIR=
export TRIM_GALORE=$SOFTWARE_DIR/TrimGalore-0.6.5/trim_galore 
export BOWTIE_TOOL=$SOFTWARE_DIR/bowtie2/bowtie2-2.0.0-beta7/bowtie2
export PHLAT_TOOL=$SOFTWARE_DIR/phlat-release
export OPTITYPE_TOOL=$SOFTWARE_DIR/OptiType/OptiTypePipeline.py
export HLAMINER_TOOL=$SOFTWARE_DIR/HLAminer/HLAminer-1.4/bin/HLAminer.pl

### Location of conda.sh
export CONDA_SH=
### See yml files in `./sequencing_pipeline/` for conda environment details

### Location of grch38 index for HISAT
export INDEX_DIR=
### Location of HLA index for HLAMINER
export HLAMINER_REF=

export N_CORES=$SLURM_CPUS_PER_TASK

# Create dirs
if [ ! -d $LOG_DIR ]; then mkdir -p $LOG_DIR;fi
if [ ! -d $TEMP_DIR ]; then mkdir -p $TEMP_DIR;fi

# SET WHICH COMPONENTS TO RUN. 0=not run, 1=run
export RUN_COMPILE=0
export RUN_PRE_QC=0
export RUN_TRIM=0
export RUN_POST_QC=0
export RUN_HISAT=0
export RUN_ARCAS=0
export RUN_PHLAT=0
export RUN_OPTITYPE=0
export RUN_HLAMINER=0

# COMPILE FASTQ FROM DIFFERENT POOLS
### NOTE: For manuscript, experiment distributed individual samples across
### multiple sequencing pools, so these sequences needed to be concatenated into
### single files for each sample, which is what this step accomplishes. For data
### already prepared as single sample files, this can be omitted, but the 
### directory names will need to be adjusted in FASTQC() and TRIM()
COMPILE_FASTQ(){
  printf "\n\
############################# COMPILE FASTQ ####################################\
  \n" >> $LOG_DIR/${1}/${1}_pipeline.log
  printf "\n### [START___FASTQ-COMPILE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  cat $SAMPLE_DIR/${1}*R2.fastq.gz > $TEMP_DIR/${1}_R2.fastq.gz
  printf "### [COMPLETE___FASTQ-COMPILE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log 
}
export -f COMPILE_FASTQ



# RUN FASTQC
# 1st argument is the sample
# 2nd argument is pre/post [0/1] trimgalore
FASTQC(){
  # This segments detects if FASTQC is being run pre- or post- trimming (based on $2)
  # It then specifies string and filename formats accordingly
  unset PHASE
  if [ $2 == 0 ]; then PHASE="PRE"; fi
  if [ $2 == 1 ]; then PHASE="POST"; fi
  FASTQC_LOG_DIR=$FASTQC_DIR/$PHASE
  if [ $2 == 0 ]; then QC_INPUT=$TEMP_DIR/${1}_R2.fastq.gz; fi
  if [ $2 == 1 ]; then QC_INPUT=$TEMP_DIR/${1}_R2_trimmed.fq.gz; fi
  [ ! -d $FASTQC_LOG_DIR ] && mkdir -p $FASTQC_LOG_DIR
  
  # This segment runs fastqc
  printf "\n\
############################# $PHASE-TRIM FASTQC ####################################\
  \n" >> $LOG_DIR/${1}/${1}_pipeline.log
  printf "\n### [START___$PHASE-FASTQC___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  fastqc -t $N_CORES -o $FASTQC_LOG_DIR $QC_INPUT \
  2> $LOG_DIR/${1}/${1}_${PHASE}QC.log
  printf "### [COMPLETE___$PHASE-FASTQC___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log 
}
export -f FASTQC

# RUN TRIM GALORE
TRIM(){
  echo "\
############################## TRIM GALORE ####################################\
  " >> $LOG_DIR/${1}/${1}_pipeline.log
  printf "\n### [START___TRIM___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  
  $TRIM_GALORE -j $N_CORES \
  -q 20 \
  --stringency 1 \
  --gzip \
  --length 20 \
  $TEMP_DIR/${1}_R2.fastq.gz \
  --output_dir $TEMP_DIR  \
  2> $LOG_DIR/${1}/${1}_trim.log
  
  printf "### [COMPLETE___TRIM___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
}
export -f TRIM


# DEFINE HISAT MAPPING PIPELINE
HISAT(){
  source $CONDA_SH
  conda activate samtools
  # This section needed to specify filename format based on trimming. Also create bam folder
  unset HISAT_INPUT
  if [ $RUN_TRIM == 0 ]; then HISAT_INPUT=$TEMP_DIR/${1}_R2.fastq.gz; fi
  if [ $RUN_TRIM == 1 ]; then HISAT_INPUT=$TEMP_DIR/${1}_R2_trimmed.fq.gz; fi
  [ ! -d $BAM_DIR ] && mkdir -p $BAM_DIR
  
  # This section runs hisat and index
  printf "\n\
########################## HISAT2 ALIGNMENT ####################################\
  \n" >> $LOG_DIR/${1}/${1}_pipeline.log
  printf "\n### [START___HISAT-AND-SORT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log

  hisat2 --dta -x $INDEX_DIR/genome \
  -U $HISAT_INPUT \
  2> $LOG_DIR/${1}/${1}_hisat.log |
  samtools sort -@ $N_CORES -O BAM > $BAM_DIR/${1}.bam

  printf "### [COMPLETE___HISAT-AND-SORT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  printf "\n### [START___INDEX___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  
  samtools index -@ $N_CORES $BAM_DIR/${1}.bam $BAM_DIR/${1}.bam.bai # Index .bam in STDIN and write to index to .bai

  printf "### [COMPLETE___INDEX___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  
  rm $TEMP_DIR/${1}*
  conda deactivate
}
export -f HISAT


# DEFINE ARCAS GENOTYPING PIPELINE
ARCAS(){
  source $CONDA_SH
  conda activate samtools
  printf "\n\
########################## ARCAS HLA TYPING ####################################\
  \n" >> $LOG_DIR/${1}/${1}_pipeline.log
  
  [ ! -d $ARCAS_DIR ] && mkdir -p $ARCAS_DIR
  
  printf "\n### [START___ARCAS-EXTRACT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  arcasHLA extract $BAM_DIR/${1}.bam -o $ARCAS_DIR/ --unmapped -t $N_CORES -v\
  2> $LOG_DIR/${1}/${1}_arcas_extract.log

  printf "### [COMPLETE___ARCAS-EXTRACT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  printf "\n### [START___ARCAS-GENOTYPE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  arcasHLA genotype $ARCAS_DIR/${1}.extracted.fq.gz -o $ARCAS_DIR/ -t $N_CORES -v\
  2> $LOG_DIR/${1}/${1}_arcas_genotype.log
  
  printf "### [COMPLETE___ARCAS-GENOTYPE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  
  rm $BAM_DIR/${1}*
  conda deactivate
}
export -f ARCAS


# DEFINE PHLAT GENOTYPING PIPELINE
PHLAT(){
  source $CONDA_SH
  conda activate phlat
  printf "\n\
########################## PHLAT TYPING ####################################\
  \n" >> $LOG_DIR/${1}/${1}_pipeline.log
  [ ! -d $PHLAT_DIR ] && mkdir -p $PHLAT_DIR
  printf "\n### [START___PHLAT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log

  python -O $PHLAT_TOOL/dist/PHLAT.py \
  -1 $ARCAS_DIR/${1}.extracted.fq.gz \
  -index $PHLAT_TOOL/b2folder \
  -b2url $BOWTIE_TOOL \
  -tag $1 \
  -e $PHLAT_TOOL \
  -o $PHLAT_DIR \
  -p $N_CORES \
  -pe 0 \
  -tmp 0 \
  1> $LOG_DIR/${1}/${1}_phlat.log 2>> $LOG_DIR/${1}/${1}_phlat.log
  
  printf "### [COMPLETE___PHLAT___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  conda deactivate
}
export -f PHLAT


# DEFINE OPTITYPE GENOTYPING PIPELINE
OPTITYPE(){
  source $CONDA_SH
  conda activate optitype
  printf "\n\
########################## OPTITYPE TYPING ####################################\
  \n" >> $LOG_DIR/${1}/${1}_pipeline.log
  [ ! -d $OPTITYPE_DIR ] && mkdir -p $OPTITYPE_DIR
  printf "\n### [START___OPTITYPE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log

  OptiTypePipeline.py -i $ARCAS_DIR/${1}.extracted.fq.gz \
  -o $OPTITYPE_DIR \
  -p $1 \
  --rna \
  -v \
  1> $LOG_DIR/${1}/${1}_optitype.log 2>> $LOG_DIR/${1}/${1}_optitype.log
  
  printf "### [COMPLETE___OPTITYPE___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  # rm $OPTITYPE_DIR/${1}_1.bam
  conda deactivate
}
export -f OPTITYPE


# DEFINE HLAMINER GENOTYPING PIPELINE
HLAMINER(){
  source $CONDA_SH
  conda activate hlaminer
  printf "\n\
########################## HLA_MINER TYPING ####################################\
  \n" >> $LOG_DIR/${1}/${1}_pipeline.log
  [ ! -d $HLAMINER_DIR ] && mkdir -p $HLAMINER_DIR
  printf "\n### [START___HLAMINER___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log

  ### Run bwa short read aligner
  bwa aln \
  -e 0 \
  -o 0 \
  $HLAMINER_REF/HLA-I_II_CDS.fasta \
  $ARCAS_DIR/${1}.extracted.fq.gz > $HLAMINER_DIR/${1}.sai \
  2> $LOG_DIR/${1}/${1}_hlaminer.log
  
  bwa samse \
  $HLAMINER_REF/HLA-I_II_CDS.fasta \
  $HLAMINER_DIR/${1}.sai \
  $ARCAS_DIR/${1}.extracted.fq.gz > $HLAMINER_DIR/${1}.sam \
  2>> $LOG_DIR/${1}/${1}_hlaminer.log
  
  ### Predict HLA
  $HLAMINER_TOOL \
  -a $HLAMINER_DIR/${1}.sam \
  -e 1 \
  -h $HLAMINER_REF/HLA-I_II_CDS.fasta \
  -p $HLAMINER_REF/hla_nom_p.txt \
  -l ${1} \
  -s 500 \
  1>> $LOG_DIR/${1}/${1}_hlaminer.log 2>> $LOG_DIR/${1}/${1}_hlaminer.log
  
  ### Clean up
  mv ./HLAminer_HPRA_${1}.csv $HLAMINER_DIR
  rm ./HLAminer_HPRA_${1}.log
  rm $HLAMINER_DIR/${1}*

  printf "### [COMPLETE___HLAMINER___$(date +'%D %X')]\n" >> $LOG_DIR/${1}/${1}_pipeline.log
  conda deactivate
}
export -f HLAMINER


# DEFINE HEADER FOR MAIN LOG FILE
LOG_HEADER(){
cat <<EOF
################################################################################
# Mapping and HLA genotyping of $1 scRNAseq
###############################################################################
# FASTQ location: $SAMPLE_DIR
# Index location: $INDEX_DIR
# Run fastq compilation status: $RUN_COMPILE
# Run pre-trim FASTQC status: $RUN_PRE_QC
# Run TRIM GALORE: $RUN_TRIM
# Run post-trim FASTQC status: $RUN_POST_QC
# Run HISAT status: $RUN_HISAT
# Run ARCAS status: $RUN_ARCAS
# Run PHLAT status: $RUN_PHLAT
# Run OPTITYPE status: $RUN_OPTITYPE
# Run HLAMINER status: $RUN_HLAMINER
EOF
}
export -f LOG_HEADER

# DEFINE PIPELINE CONTROLLER FUNCTION
PIPELINE(){
  [ ! -d $LOG_DIR/$1 ] && mkdir -p $LOG_DIR/$1
  echo "START: sample $1 at $(date +'%D %X')"
  LOG_HEADER $1 > $LOG_DIR/${1}/${1}_pipeline.log # CREATE LOG FILE WITH CUSTOM HEADER (SEE BELOW)
  [ $RUN_COMPILE == 1 ] && COMPILE_FASTQ $1
  [ $RUN_PRE_QC == 1 ] && FASTQC $1 0
  [ $RUN_TRIM == 1 ] && TRIM $1 
  [ $RUN_POST_QC == 1 ] && FASTQC $1 1
  [ $RUN_HISAT == 1 ] && HISAT $1
  [ $RUN_ARCAS == 1 ] && ARCAS $1
  [ $RUN_PHLAT == 1 ] && PHLAT $1
  [ $RUN_OPTITYPE == 1 ] && OPTITYPE $1
  [ $RUN_HLAMINER == 1 ] && HLAMINER $1
  echo "COMPLETE: sample $1 at $(date +'%D %X')"
}
export -f PIPELINE

# ACTUAL RUN COMMAND
cat $SAMPLE_LIST | parallel --delay 15 -j $SLURM_NTASKS --joblog $LOG_DIR/parallel.log PIPELINE {}
