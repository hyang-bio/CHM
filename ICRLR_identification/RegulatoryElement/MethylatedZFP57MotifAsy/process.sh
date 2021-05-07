#!/bin/bash

# ———————————————————————————————————
# Obtain score of allele-specific methylated ZFP57 motif in asCHM.
# ———————————————————————————————————

# Obtain location of ZFP57 motif in asCHM
python getSeqLoc.py "TGCCGC" asCHM.bed asCHM.ZFP57motif.bed mm10.2bit

# DNA methylation level in covered ZFP57 motif
echo -ne "#Chrom\tStart\tEnd" > title.txt
cut -f 1-3 asCHM.ZFP57motif.bed | sort -k1,1 -k2,2n > asCHM.ZFP57motif.methyl_merged.txt
for sample in Oocyte Sperm 2cell.GG 2cell.AG Morula.GG Morula.AG ICM.GG ICM.AG TE.GG TE.AG
do
	bash averageMethylInRegionMultipleThreads.sh \
	 asCHM.ZFP57motif.bed \
	 ${sample}.sam.G.bed \
	 asCHM.ZFP57motif.methyl_${sample}.txt \
	 10

	echo -ne "\t${sample}" >> title.txt
	cut -f 4 asCHM.ZFP57motif.methyl_${sample}.txt | paste asCHM.ZFP57motif.methyl_merged.txt - > tmp && mv tmp asCHM.ZFP57motif.methyl_merged.txt
done # for sample end

# Count samples with asymmetric methylated ZFP57 motif
cat asCHM.ZFP57motif.methyl_merged.txt | awk 'function score_M(M_mat, M_pat){\
 if(M_mat!="NA" && M_pat!="NA" && M_mat>=0.5 && M_pat/M_mat<=0.5){return 1} else if(M_mat!="NA" && M_pat!="NA" && M_pat>=0.5 && M_mat/M_pat<=0.5){return -1} else{return 0}};\
 BEGIN{FS=OFS="\t"}{print $0, score_M($4, $5)+score_M($6, $7)+score_M($8, $9)+score_M($10, $11)+score_M($12, $13)}' > asCHM.ZFP57motif.methyl_merged.score.txt