#!/bin/bash
#slurm options
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH -q normal
#SBATCH -J ccFDR_RA_AF
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH -o ccFDR_RA_AF.log

cd /storage/zhenghoufengLab/qianyu/Software/pleioFDR/pleiofdr

module load matlab/R2021a
matlab -nodisplay -nosplash < runme.m
