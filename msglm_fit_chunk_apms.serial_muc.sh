#!/bin/bash
#SBATCH -o /gpfs/scratch/pn69ha/ge68wan2/logs/cov2/msglm_fit_%A_%a.%N.log
#SBATCH -J ast_cov2_msglm_apms
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --mincpus=8
##SBATCH --mem-per-cpu=2GB
##SBATCH --exclusive=user
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-socket=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexey.stukalov@tum.de
#SBATCH --export=NONE
#SBATCH --time=72:00:00

source /etc/profile.d/modules.sh

#module switch spack/release/19.1
module load charliecloud
#export PATH=/lrz/sys/tools/charliecloud/0.9.9/bin:$PATH

PROJECT_ID=cov2
DATA_VERSION=20200329
FIT_VERSION=20200329
CHUNK_IDS_FILE=${SCRATCH}/${PROJECT_ID}/pending_chunk_ids
if [[ -f $CHUNK_IDS_FILE ]]; then
  echo "Reading ${CHUNK_IDS_FILE}..."
  readarray -t CHUNK_IDS < $CHUNK_IDS_FILE
fi

CHUNK_IX_START=$((4000+SLURM_ARRAY_TASK_ID))
CHUNK_IX_LAST=5542
# chunk increment, should be equal to the number of jobs in array
CHUNK_IX_STRIDE=$SLURM_ARRAY_TASK_COUNT
#CHUNK_IX_STRIDE=$SLURM_ARRAY_TASK_MAX
# if CHUNK_IDS are set, use their length
if [[ -v CHUNK_IDS ]]; then
  CHUNK_IX_LAST=${#CHUNK_IDS[@]}
fi
echo "Chunks indices to process: ${CHUNK_IX_START}..${CHUNK_IX_LAST}..${CHUNK_IX_STRIDE}"

# chunks farming, sequentially run CHUNK_IX_LAST/CHUNK_IX_STRIDE jobs
for ((chunk_ix = CHUNK_IX_START; chunk_ix <= CHUNK_IX_LAST; chunk_ix += CHUNK_IX_STRIDE)); do
if [[ -v CHUNK_IDS ]]; then
  chunk_id=${CHUNK_IDS[((chunk_ix-1))]}
else
  chunk_id=$chunk_ix
fi
echo "Chunk #$chunk_id: MSGLM fit starting..."
{
ch-run $SCRATCH/docker4muc/archpc.msglm \
  -t --unset-env='*PATH' \
  -b $HOME/projects/adhoc:/projects/adhoc \
  -b $HOME/data:/data \
  -b $HOME/analysis:/analysis \
  -b $SCRATCH:/scratch \
  -- Rscript /projects/adhoc/$PROJECT_ID/msglm_fit_chunk_apms.R \
  $PROJECT_ID $SLURM_JOB_NAME mq_apms_20200329 $DATA_VERSION $FIT_VERSION $SLURM_ARRAY_JOB_ID $chunk_id && \
echo "Chunk #$chunk_id: MSGLM fit done"
} || {
echo "Chunk #$chunk_id: MSGLM fit failed"
}
done
echo "All chunks of array task id=$SLURM_ARRAY_TASK_ID finished"

