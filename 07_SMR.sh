#!/bin/bash
#slurm options
#SBATCH -p amd-ep2,intel-e5
#SBATCH -q normal
#SBATCH -J AF_smr_v8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=20G
#SBATCH -o AF_smr_v8.log

cd /storage/zhenghoufengLab/qianyu/project/RA_CVD/SMR/

/storage/zhenghoufengLab/qianyu/Software/SMR/smr_Linux \
--bfile /storage/zhenghoufengLab/qianyu/Software/g1000_eur \
--gwas-summary /storage/zhenghoufengLab/qianyu/project/RA_CVD/SMR/AF_SMR.ma \
--beqtl-summary /storage/zhenghoufengLab/share/database/eqtl/GTEx_V8_cis_eqtl_summary_hg19/Artery_Aorta/Artery_Aorta \
--peqtl-smr 5e-8 \
--ld-upper-limit 0.9 \
--ld-lower-limit 0.05 \
--peqtl-heidi 1.57e-3 \
--heidi-min-m 3 \
--heidi-max-m 20 \
--cis-wind 2000 \
--thread-num 5 \
--out /storage/zhenghoufengLab/qianyu/project/RA_CVD/SMR/AF_Artery_Aorta
	
