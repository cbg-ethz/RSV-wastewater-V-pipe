#!/bin/bash

#SBATCH --job-name=lollipop_rsvb       # Job name
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem-per-cpu=2G            # Memory per node
#SBATCH --time=72:00:00             # Maximum runtime
module load eth_proxy

# Reset SECONDS to 0
SECONDS=0

ldata="./rsv_b_lofreq_2022_2023"

lollipop deconvolute $ldata/tallymut_EPI_ISL_1653999.tsv \
    -o $ldata/deconvolved.csv \
    --variants-config $ldata/variant_config.yaml \
    --variants-dates $ldata/var_dates.yaml \
    --deconv-config $ldata/deconv_config.yaml \
    --seed=42 \
    --n-cores=8 \
    > $ldata/deconvolution.log 2>&1

# Print elapsed time
echo "Time taken: $SECONDS seconds"
