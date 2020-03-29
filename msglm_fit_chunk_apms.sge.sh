#!/bin/bash
#$ -pe openmp 8
#$ -R y
#$ -l h_rt=72:00:00
#$ -l h_vmem=12G
#$ -N ast_cov2_apms_msglm
##$ -S /fs/home/stukalov/gentoo/bin/bash
#$ -m ae
#$ -j y
##$ -w v
#$ -o /fs/pool/pool-innate-analysis/scratch/stukalov/logs/cov2/$JOB_NAME_$JOB_ID_$TASK_ID.log

IMAGES_PATH=/fs/pool/pool-innate-bioinfo2/ccimages
FS_PREFIX=/fs/pool/pool-innate-analysis/stukalov
PROJECT_ID=cov2

# disable MKL threading
export MKL_NUM_THREADS=1
#export TMPDIR=$HOME/scratch/Rtmp_qsub
#JOB_ID=1
#JOB_NAME="ast_cov2_apms_msglm"
#SGE_TASK_ID=3
#export NSLOTS=8

module load charliecloud

# do it!
ch-run $IMAGES_PATH/archpc.msglm \
  -t --unset-env='*PATH' \
  -b $HOME/projects/adhoc:/projects/adhoc \
  -b $FS_PREFIX/analysis:/analysis \
  -b $FS_PREFIX/scratch:/scratch \
  -- Rscript /projects/adhoc/$PROJECT_ID/msglm_fit_chunk_apms.R \
  $PROJECT_ID $JOB_NAME mq_apms_20200329 20200329 20200329 $JOB_ID $SGE_TASK_ID
