source activate ldsc13

cd /storage/zhenghoufengLab/qianyu/data/Sumstats/Brain_and_heart_disease

/storage/zhenghoufengLab/qianyu/Software/ldsc/ldsc.py \
--rg RA.sumstats.gz,HF.sumstats.gz \
--samp-prev 0.246,0.048 \
--pop-prev 0.019,0.0202 \
--ref-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--w-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--out /storage/zhenghoufengLab/qianyu/project/RA_CVD/gen_corr/RA_HF_h2

/storage/zhenghoufengLab/qianyu/Software/ldsc/ldsc.py \
--rg RA.sumstats.gz,AF.sumstats.gz \
--samp-prev 0.246,0.058806638 \
--pop-prev 0.019,0.04972315 \
--ref-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--w-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--out /storage/zhenghoufengLab/qianyu/project/RA_CVD/gen_corr/RA_AF_h2

/storage/zhenghoufengLab/qianyu/Software/ldsc/ldsc.py \
--rg RA.sumstats.gz,CAD.sumstats.gz \
--samp-prev 0.246,0.223046429 \
--pop-prev 0.019,0.08234054 \
--ref-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--w-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--out /storage/zhenghoufengLab/qianyu/project/RA_CVD/gen_corr/RA_CAD_h2

/storage/zhenghoufengLab/qianyu/Software/ldsc/ldsc.py \
--rg RA.sumstats.gz,MI_2.sumstats.gz \
--samp-prev 0.246,0.254117866 \
--pop-prev 0.019,0.04454882 \
--ref-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--w-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--out /storage/zhenghoufengLab/qianyu/project/RA_CVD/gen_corr/RA_MI_h2

/storage/zhenghoufengLab/qianyu/Software/ldsc/ldsc.py \
--rg RA.sumstats.gz,AS.sumstats.gz \
--samp-prev 0.246,0.090855974 \
--pop-prev 0.019,0.02906972 \
--ref-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--w-ld-chr /storage/zhenghoufengLab/qianyu/Software/ldsc/eur_w_ld_chr/ \
--out /storage/zhenghoufengLab/qianyu/project/RA_CVD/gen_corr/RA_stroke_h2
