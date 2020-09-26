#! /bin/sh
#SBATCH -o /gpfs/scratch/pn69ha/ge68wan2/logs/cov2/msglm_fit_parsars_ptm_%A_%a.log
#SBATCH -J cov2_msglm_parsars_ptm
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_large
#SBATCH --reservation=special_covid
#SBATCH --nodes=60
##SBATCH --ntasks=1
#SBATCH --mincpus=7
##SBATCH --mem-per-cpu=2GB
##SBATCH --exclusive=user
#SBATCH --cpus-per-task=7
##SBATCH --ntasks-per-socket=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexey.stukalov@tum.de
#SBATCH --export=NONE
#SBATCH --time=192:00:00

source /etc/profile.d/modules.sh
module load slurm_setup
module load charliecloud

PROJECT_ID=cov2
DATA_VERSION=20200920
FIT_VERSION=20200920
MS_FOLDER=snaut_parsars_ptm_20200907
CHUNKS_PER_JOB=10000
MAX_CHUNK=27240
CHUNK_IDS_FILE=${SCRATCH}/${PROJECT_ID}/${PROJECT_ID}_${MS_FOLDER}_${FIT_VERSION}_pending_chunk_ids
if [[ -f $CHUNK_IDS_FILE ]]; then
  echo "Reading ${CHUNK_IDS_FILE}..."
  readarray -t CHUNK_IDS < $CHUNK_IDS_FILE
else
  CHUNK_IDS=( $(seq 1 $MAX_CHUNK) )
fi

CHUNK_IDS=(${CHUNK_IDS[@]:$(((SLURM_ARRAY_TASK_ID-1)*CHUNKS_PER_JOB)):$CHUNKS_PER_JOB})

for chunk in ${CHUNK_IDS[@]}
do
srun -c7 -n1 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=2GB --wait=0 \
     --no-kill --distribution=block --exclusive=user \
ch-run $SCRATCH/docker4muc/archpc.msglm \
  -t --unset-env='*PATH' \
  -b $HOME/projects/adhoc:/projects/adhoc \
  -b $HOME/data:/data \
  -b $HOME/analysis:/analysis \
  -b $SCRATCH:/scratch \
  -- Rscript /projects/adhoc/$PROJECT_ID/msglm_fit_chunk_parsars_ptm.R \
  $PROJECT_ID $SLURM_JOB_NAME $MS_FOLDER $DATA_VERSION $FIT_VERSION $SLURM_JOB_ID $chunk &
# sleep to avoid clogging of too many steps
if ! ((chunk % 300)); then
  sleep 1
fi
done

wait
echo "All chunks of job #$SLURM_JOB_ID finished"

