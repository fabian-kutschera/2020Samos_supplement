#!/bin/bash
#SBATCH -J Samos_Ryo_noWL_v3
#SBATCH -o /hppfs/scratch/09/ru64lev2/Samos_final_Ryo_v3/%x.%j.out
#SBATCH -e /hppfs/scratch/09/ru64lev2/Samos_final_Ryo_v3/%x.%j.out

#Initial working directory:
#SBATCH --chdir=./

#Notification and type
#SBATCH --mail-type=BEGIN,END,ARRAY_TASKS
#SBATCH --mail-user=f.kutschera@campus.lmu.de

#SBATCH --ear=off
#SBATCH --time=14:30:00
#SBATCH --no-requeue
#SBATCH --export=ALL
#SBATCH --account=pn49ha
#SBATCH --partition=general
#SBATCH --nodes=256
#SBATCH --ntasks-per-node=1
module load slurm_setup

#Run the program:
export MP_SINGLE_THREAD=no
unset KMP_AFFINITY
export OMP_NUM_THREADS=94
export OMP_PLACES="cores(47)"
#Prevents errors such as experience in Issue #691
export I_MPI_SHM_HEAP_VSIZE=8192

export XDMFWRITER_ALIGNMENT=8388608
export XDMFWRITER_BLOCK_SIZE=8388608
export SC_CHECKPOINT_ALIGNMENT=8388608

export SEISSOL_CHECKPOINT_ALIGNMENT=8388608
export SEISSOL_CHECKPOINT_DIRECT=1
export ASYNC_MODE=THREAD
export ASYNC_BUFFER_ALIGNMENT=8388608
source /etc/profile.d/modules.sh

echo 'num_nodes:' $SLURM_JOB_NUM_NODES 'ntasks:' $SLURM_NTASKS 'cpus_per_task:' $SLURM_CPUS_PER_TASK
ulimit -Ss 2097152

mpiexec -n $SLURM_NTASKS /hppfs/work/pn49ha/ru64lev2/SeisSol/build-release/SeisSol_Release_dskx_5_elastic parameters_samos_noWL.par
