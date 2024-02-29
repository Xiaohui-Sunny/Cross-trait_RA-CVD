cat /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/pheno.txt |\
	while read i
	do
		for m in $(seq 1 20)
		do
		echo "#!/bin/bash
#slurm options
#SBATCH -p intel-sc3
#SBATCH -q huge
#SBATCH -J mixer_S1_${i}_${m}
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH -o mixer_S1_${i}_${m}.log

source activate mixer

python /home/zhenghoufengLab/qianyu/mixer/precimed/mixer.py fit1 \\
      --trait1-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/${i}_mixer.txt.gz \\
      --extract /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${m}.snps \\
      --bim-file /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \\
      --ld-file /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \\
      --lib  /home/zhenghoufengLab/qianyu/mixer/src/build/lib/libbgmg.so \\
      --out /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/${i}.fit.rep${m}
" > /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S1_${i}_${m}.sh
sbatch /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S1_${i}_${m}.sh
done
done


cat  /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/trait.txt|\
	while read i
	do
		for m in $(seq 1 20)
		do
		echo "#!/bin/bash
#slurm options
#SBATCH -p intel-sc3
#SBATCH -q huge
#SBATCH -J mixer_S2_${i}_${m}
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH -o mixer_S2_${i}_${m}.log

source activate mixer

python /home/zhenghoufengLab/qianyu/mixer/precimed/mixer.py fit2 \\
      --trait1-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/RA_mixer.txt.gz \\
      --trait2-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/${i}_mixer.txt.gz \\
      --trait1-params-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA.fit.rep${m}.json \\
      --trait2-params-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/${i}.fit.rep${m}.json \\
      --extract /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${m}.snps \\
      --bim-file /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \\
      --ld-file /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \\
      --lib  /home/zhenghoufengLab/qianyu/mixer/src/build/lib/libbgmg.so \\
      --out /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.fit.rep${m}
" > /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S2_RA_${i}_${m}.sh
sbatch /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S2_RA_${i}_${m}.sh
done
done

cat /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/trait.txt |\
	while read i
	do
		for m in $(seq 1 20)
		do
		echo "#!/bin/bash
#slurm options
#SBATCH -p intel-sc3
#SBATCH -q huge
#SBATCH -J mixer_S3_${i}_${m}
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH -o mixer_S3_${i}_${m}.log

source activate mixer

python /home/zhenghoufengLab/qianyu/mixer/precimed/mixer.py test2 \\
      --trait1-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/RA_mixer.txt.gz \\
      --trait2-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/${i}_mixer.txt.gz \\
      --load-params-file /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.fit.rep${m}.json \\
      --bim-file /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \\
      --ld-file /storage/zhenghoufengLab/niuyuxiao/software/MiXeR/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \\
      --lib  /home/zhenghoufengLab/qianyu/mixer/src/build/lib/libbgmg.so \\
      --out /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.test.rep${m}
" > /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S3_RA_${i}_${m}.sh
sbatch /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S3_RA_${i}_${m}.sh
done
done


cat /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/trait.txt |\
	while read i
	do
		echo "#!/bin/bash
#slurm options
#SBATCH -p intel-sc3
#SBATCH -q huge
#SBATCH -J mixer_S4_${i}
#SBATCH -c 6
#SBATCH --mem=10G
#SBATCH -o mixer_S4_${i}.log

source activate mixer

python /home/zhenghoufengLab/qianyu/mixer/precimed/mixer_figures.py combine \\
      --json /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.fit.rep@.json \\
      --out /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.fit

python /home/zhenghoufengLab/qianyu/mixer/precimed/mixer_figures.py combine \\
      --json /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.test.rep@.json \\
      --out /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.test

python /home/zhenghoufengLab/qianyu/mixer/precimed/mixer_figures.py two \\
      --json-fit /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.fit.json \\
      --json-test /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}.test.json \\
      --statistic mean std \\
      --out /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/RA_vs_${i}
" > /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S4_RA_${i}.sh
sbatch /storage/zhenghoufengLab/qianyu/project/RA_CVD/MixeR/data/shell/mixer_S4_RA_${i}.sh
done
