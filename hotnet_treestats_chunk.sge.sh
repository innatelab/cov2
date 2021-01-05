#!/bin/bash
#$ -pe openmp 2
#$ -R y
#$ -l h_rt=48:00:00
#$ -l h_vmem=32G
#$ -N cov2_hotnet_real
#$ -m ae
#$ -j y
##$ -w v
#$ -o /fs/pool/pool-innate-analysis/scratch/stukalov/logs/cov2/$JOB_NAME_$JOB_ID_$TASK_ID.log

#IMAGES_PATH=/fs/pool/pool-innate-bioinfo2/ccimages
IMAGES_PATH=/fs/home/stukalov/tmp_ccimages
FS_PREFIX=/fs/pool/pool-innate-analysis/stukalov
CC_PATH=/fs/home/stukalov/apps/charliecloud/bin

PROJECT_ID=cov2
HOTNET_VERSION=20201022

chunk=${SGE_TASK_ID}

#module load charliecloud

echo "Starting the analysis of chunk #${chunk}..."

export MKL_NUM_THREADS=$NSLOTS

# do it!
$CC_PATH/ch-run $IMAGES_PATH/archpc.julia \
     -t --no-home \
     --unset-env=LD_LIBRARY_PATH --set-env=$HOME/tmp_projects/adhoc/$PROJECT_ID/hotnet_sge.julia.env \
     -b $HOME/tmp_projects/adhoc:/projects/adhoc \
     -b $FS_PREFIX/scratch:/scratch -- \
     julia --project=/projects/adhoc/$PROJECT_ID \
     /projects/adhoc/$PROJECT_ID/hotnet_treestats_chunk.jl \
     $PROJECT_ID $JOB_NAME $HOTNET_VERSION $JOB_ID $chunk 10

echo "Chunk analysis finished"
