#!/bin/bash
#SBATCH --job-name=TeloNPRun
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --partition=cpu2023
#SBATCH --mem=250GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ramin.kahidi@ucalgary.ca

source ~/software/init-conda
conda activate pytorch

cd "/home/ramin.kahidi/GreiderLab/TeloNP/TeloBP/Demo/"

echo ==============================================
echo Hostname: `hostname`
echo Job ID: $SLURM_JOBID

echo Node: $SLURM_NODELIST

start=`date +%Y-%m-%d\ %H:%M:%S`

python TeloNPScript.py

end=`date +%Y-%m-%d\ %H:%M:%S`

echo ==============================================
echo Start: $start
echo End: $end
