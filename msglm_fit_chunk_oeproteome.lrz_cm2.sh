#! /bin/sh
#SBATCH -o /gpfs/scratch/pn69ha/ge68wan2/logs/cov2/msglm_fit_oeprot_%J.log
#SBATCH -J ast_cov2_msglm_apms
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_large
#SBATCH --reservation=special_covid
#SBATCH --nodes=108
##SBATCH --ntasks=1
#SBATCH --mincpus=7
##SBATCH --mem-per-cpu=2GB
##SBATCH --exclusive=user
#SBATCH --cpus-per-task=7
##SBATCH --ntasks-per-socket=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexey.stukalov@tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

source /etc/profile.d/modules.sh
module load slurm_setup
module load charliecloud

PROJECT_ID=cov2
DATA_VERSION=20200411
FIT_VERSION=20200411
MS_FOLDER=spectronaut_oeproteome_20200411
CHUNK_IDS_FILE=${SCRATCH}/${PROJECT_ID}/${PROJECT_ID}_${MS_FOLDER}_${FIT_VERSION}_pending_chunk_ids
if [[ -f $CHUNK_IDS_FILE ]]; then
  echo "Reading ${CHUNK_IDS_FILE}..."
  readarray -t CHUNK_IDS < $CHUNK_IDS_FILE
else
  CHUNK_IDS=( $(seq 1 6907) )
fi

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
  -- Rscript /projects/adhoc/$PROJECT_ID/msglm_fit_chunk_oeproteome.R \
  $PROJECT_ID $SLURM_JOB_NAME $MQ_FOLDER $DATA_VERSION $FIT_VERSION $SLURM_JOB_ID $chunk &
# sleep to avoid clogging of too many steps
if ! ((chunk % 400)); then
  sleep 1
fi
done

wait
echo "All chunks of job #$SLURM_JOB_ID finished"

