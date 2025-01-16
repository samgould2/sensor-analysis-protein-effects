#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of tasks (really number of CPU Cores/task). Don't specify more than 16 unless approved by the system admin
#SBATCH --array=1-60    #change this to match up with the number of parallel computing jobs; also generate a config file 
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=samgould@mit.edu           # Email to which notifications will be sent. Equivalent to the -M option in SGE. You must replace [] with your email address.
#SBATCH --exclude=c[5-22]
#SBATCH --nice=100000
#############################################

module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit

#create a new virtual environment if needed; this one is fine from previously
conda activate /home/samgould/.conda/envs/sensor_lib_sg

#go to correct directory
cd /net/bmc-lab2/data/lab/sanchezrivera/samgould/241105SanA

#access the config file
#NEED TO CHANGE THE CONFIG FILE FOR YOUR SPECIFIC RUN
config=./CDK_config.txt

# Extract R1_FILE name for the current $SLURM_ARRAY_TASK_ID
R1_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
R2_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
folder_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

EDITOR='CBE'

python3 sensor_HGVSp_pipeline_v2.py master_df_HGVSp_quant_pool1_v2.csv cdks_info_updated.csv ${EDITOR} ${folder_name}
