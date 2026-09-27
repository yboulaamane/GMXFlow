#!/bin/bash
# Example SLURM job for the production run. The #SBATCH lines are for one
# particular cluster: change the account, QOS and partition to yours.
#SBATCH --account=gpu_users
#SBATCH --qos=gpu
#SBATCH --job-name=md_run            # Job name
#SBATCH --ntasks=16                  # Total number of tasks (CPUs)
#SBATCH --cpus-per-task=1            # CPUs per task
#SBATCH --gres=gpu:1                 # Request 1 GPU
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --time=7-00:00:00            # Maximum runtime (7 days)
#SBATCH --partition=gpu-prodq        # GPU partition
#SBATCH --output=md_run_%j.out       # Standard output (%j is the job ID)
#SBATCH --error=md_run_%j.err        # Standard error (%j is the job ID)

# Load GROMACS if your cluster uses environment modules, e.g.
#   GMX_MODULE=GROMACS/2021.3-foss-2021a-CUDA-11.3.1 sbatch submit_gromacs_gpu.sh
if [ -n "${GMX_MODULE:-}" ]; then
    module load "$GMX_MODULE"
fi

# Log allocated resources
echo "Job running on node(s): $SLURM_JOB_NODELIST"
echo "Allocated GPU(s): $SLURM_GPUS"
echo "CPU(s) per task: $SLURM_CPUS_PER_TASK"

# Run GROMACS simulation 
"${GMX:-gmx}" mdrun -deffnm "${DEFFNM:-md_0_100}"
