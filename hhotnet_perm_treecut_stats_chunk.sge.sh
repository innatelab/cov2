#!/bin/bash
#$ -pe openmp 8
#$ -R y
#$ -l h_rt=48:00:00
#$ -l h_vmem=24G
#$ -N cov2_hotnet_stats
##$ -S /fs/home/stukalov/gentoo/bin/bash
#$ -m ae
#$ -j y
##$ -w v
#$ -o /fs/pool/pool-innate-analysis/scratch/stukalov/logs/cov2/$JOB_NAME_$JOB_ID_$TASK_ID.log

IMAGES_PATH=/fs/pool/pool-innate-bioinfo2/ccimages
FS_PREFIX=/fs/pool/pool-innate-analysis/stukalov
PROJECT_ID=cov2
HOTNET_VERSION=20200607

chunk=${SGE_TASK_ID}

module load charliecloud

echo "Starting the analysis of chunk #${chunk}..."

# do it!
ch-run $IMAGES_PATH/archpc.julia \
     -t --no-home \
     --unset-env=LD_LIBRARY_PATH --set-env=$HOME/projects/adhoc/$PROJECT_ID/archpc.julia.env \
     -b $HOME/projects/adhoc:/projects/adhoc \
     -b $FS_PREFIX/scratch:/scratch -- \
     julia --project=/projects/adhoc/$PROJECT_ID \
     /projects/adhoc/$PROJECT_ID/hhotnet_perm_treecut_stats_chunk.jl \
     $PROJECT_ID $JOB_NAME $HOTNET_VERSION $JOB_ID $chunk

echo "Chunk analysis finished"
