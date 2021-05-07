#!/bin/bash


# ———————————————————————————————————
# Obtain allele-specific expression pattern of gene or transposable element
# according to expression.
# ———————————————————————————————————


# Gene
cat mm10.refGene.log2Plus1_AveReps_FPKM.txt | awk 'function max(a, b){return a > b ? a : b};\
 function exprPattern(fpkm_mat, fpkm_pat){if(fpkm_mat >= 2 && fpkm_mat-fpkm_pat >= 1){return 1} else if(fpkm_pat >= 2 && fpkm_pat-fpkm_mat >= 1){return -1} else{return 0}};\
 BEGIN{FS=OFS="\t"}{\
 if(NR==1){print "#Chrom", "Start", "End", "RefSeq", "GeneSymbol", "Strand", "2-cell", "Morula", "ICM", "TE"} else{\
 	if($6=="+"){print $1, max(0, $2-2000), $2+2000, $4, $5, $6, exprPattern($7, $8), exprPattern($9,$10), exprPattern($11, $12), exprPattern($13, $14)} else{\
 	print $1, max(0, $3-2000), $3+2000, $4, $5, $6, exprPattern($7, $8), exprPattern($9,$10), exprPattern($11, $12), exprPattern($13, $14)}}}' | \
 (sed -u 1q;sort -k1,1 -k2,2n) > mm10.refGene.exprPattern.txt

 # Transposable element
 cat mm10.te_tx.log2Plus1_AveReps_FPKM.txt | awk 'function max(a, b){return a > b ? a : b};\
 function exprPattern(fpkm_mat, fpkm_pat){if(fpkm_mat >= 2 && fpkm_mat-fpkm_pat >= 1){return 1} else if(fpkm_pat >= 2 && fpkm_pat-fpkm_mat >= 1){return -1} else{return 0}};\
 BEGIN{FS=OFS="\t"}{\
 if(NR==1){print "#Chrom", "Start", "End", "TransposableElementTranscript", "TransposableElementSymbol", "Strand", "2-cell", "Morula", "ICM", "TE"} else{\
 	if($6=="+"){print $1, max(0, $2-2000), $2+2000, $4, $5, $6, exprPattern($7, $8), exprPattern($9,$10), exprPattern($11, $12), exprPattern($13, $14)} else{\
 	print $1, max(0, $3-2000), $3+2000, $4, $5, $6, exprPattern($7, $8), exprPattern($9,$10), exprPattern($11, $12), exprPattern($13, $14)}}}' | \
 (sed -u 1q;sort -k1,1 -k2,2n) > mm10.te_tx.exprPattern.txt