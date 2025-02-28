#!/bin/bash
#SBATCH --job-name=vpipe_rsv       # Job name
#SBATCH --ntasks=1                  # Number of tasks (usually 1 for V-pipe)
#SBATCH --cpus-per-task=16           # Number of CPU cores per task
#SBATCH --mem-per-cpu=8G            # Memory per node
#SBATCH --time=96:00:00             # Maximum runtime
module load eth_proxy
# Run V-pipe
./vpipe -p --cores 16 --rerun-incomplete
