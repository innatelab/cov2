#! /bin/sh
#SBATCH -J ast_cov2_hotnet_perm
#SBATCH -o /gpfs/scratch/pn69ha/ge68wan2/logs/cov2/hotnet_perm_%J.log
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_large
#SBATCH --reservation=special_covid
#SBATCH --nodes=60
##SBATCH --ntasks=1
#SBATCH --mincpus=1
##SBATCH --mem-per-cpu=18GB
##SBATCH --exclusive=user
#SBATCH --cpus-per-task=1
##SBATCH --ntasks-per-socket=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexey.stukalov@tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

source /etc/profile.d/modules.sh
module load slurm_setup
module load charliecloud

PROJECT_ID=cov2
HOTNET_VERSION=20200917

# treestats for real data
for chunk in {1..6}
do
srun -c1 -n1 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=14GB --wait=0 \
     --no-kill --distribution=cyclic --exclusive=user \
     ch-run -t --no-home \
     --unset-env=LD_LIBRARY_PATH --set-env=$HOME/projects/adhoc/$PROJECT_ID/hotnet_cm2lrz.julia.env \
     -b $HOME/projects/adhoc:/projects/adhoc \
     -b $SCRATCH:/scratch \
     $SCRATCH/docker4muc/archpc.julia -- \
     julia --project=/projects/adhoc/$PROJECT_ID \
     /projects/adhoc/$PROJECT_ID/hotnet_treestats_chunk.jl \
     $PROJECT_ID $SLURM_JOB_NAME $HOTNET_VERSION $SLURM_JOB_ID $chunk 10 &
if !((chunk % 200)); then
  sleep 1
fi
done

# permutations
for chunk in {1..1200}
do
srun -c1 -n1 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=14GB --wait=0 \
     --no-kill --distribution=cyclic --exclusive=user \
     ch-run -t --no-home \
     --unset-env=LD_LIBRARY_PATH --set-env=$HOME/projects/adhoc/$PROJECT_ID/hotnet_cm2lrz.julia.env \
     -b $HOME/projects/adhoc:/projects/adhoc \
     -b $SCRATCH:/scratch \
     $SCRATCH/docker4muc/archpc.julia -- \
     julia --project=/projects/adhoc/$PROJECT_ID \
     /projects/adhoc/$PROJECT_ID/hotnet_perm_chunk.jl \
     $PROJECT_ID $SLURM_JOB_NAME $HOTNET_VERSION $SLURM_JOB_ID $chunk 50 &
if !((chunk % 200)); then
  sleep 1
fi
done

wait
