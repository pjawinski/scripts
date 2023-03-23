#!/bin/bash

# =======================================================================
# === code snippets for plotting height vs. polygenic score of height ===
# =======================================================================
# feel free to adapt these code snippets for your own purposes
# note that this script calls pgsheight.R, which needs to be downloaded and adapted, too!

# set working directory
cd /your/working/directory/

# download PGS weights (see https://www.joelhirschhornlab.org/giant-consortium-results)
mkdir -p data/pgs/
wget -P data/pgs/ --content-disposition https://504394d8-624a-4827-9f25-95a83cd9675a.filesusr.com/archives/1d0101_d26d04152fd048ec8cf9787c6f1a59d2.gz?dn=GIANT_HEIGHT_YENGO_2022_PGS_WEIGHTS_ALL.gz

# calculate polygenic score for height
mkdir -p results/height/pgs/chr
for i in {1..22} X; do
	plink2 \
		--pfile data/genetics/processed/imp_qc/chr${i} \
		--score data/pgs/GIANT_HEIGHT_YENGO_2022_PGS_WEIGHTS_ALL.gz 1 5 6 header cols=+scoresums \
		--out results/height/pgs/chr/chr${i}
done

awk 'NR==1 { print "FID", "IID", "CHR_COUNT", "NMISS_ALLELE_CT", "SCORE1_SUM"; next }
	 NR==FNR { id[$1]=$1; chrcount[$1]=1; nmiss[$1]=$3; score[$1]=$6; next }
	 FNR==1 { next }
	 { chrcount[$1]++; nmiss[$1]=nmiss[$1]+$3; score[$1]=score[$1]+$6 }
	 END { for (i in score) { print id[i], id[i], chrcount[i], nmiss[i], score[i] } }' OFS='\t' results/height/pgs/chr/*.sscore \
	 > results/height/pgs/pgs.score

# create pgs plot
Rscript code/pgsheight.R
