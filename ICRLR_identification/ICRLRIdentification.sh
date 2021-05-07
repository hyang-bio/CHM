#!/bin/bash

# ———————————————————————————————————
# To identify ICRLR with transcriptional regulation functions similar to those of known ICR,
# we calculated scores for all identified allele-specific CHMs based on features of known ICRs.
# Three groups of features were used, including epigenetic modification asymmetry, regulatory elements, 
# and surrounding genes and transposable elements with allele-specific expression. 
# ———————————————————————————————————

dirPATH=$(pwd)
epiAsym_F(){
	mkdir -p ${dirPATH};cd ${dirPATH}
	
	# Manage DNA methylation level/H3K9me3 signal in each sample
	echo -ne "#Chrom\tStart\tEnd\tName" > title.txt
	cat asCHM.bed > asCHM.methyl_merged.txt
	cat asCHM.bed > asCHM.H3K9me3_merged.txt
	for sample in Oocyte Sperm 2cell.GG 2cell.AG Morula.GG Morula.AG ICM.GG ICM.AG TE.GG TE.AG
	do
		echo -ne "\t${sample}" >> title.txt
		cut -f 4 ${dirPATH}/Methyl/asCHM.methyl_${sample}.txt | paste asCHM.methyl_merged.txt - > tmp && mv tmp asCHM.methyl_merged.txt
		cut -f 4 ${dirPATH}/Methyl/asCHM.H3K9me3_${sample}.txt | paste asCHM.H3K9me3_merged.txt - > tmp && mv tmp asCHM.H3K9me3_merged.txt
	done # for sample end
	echo -ne "\n" >> title.txt
	cat title.txt asCHM.methyl_merged.txt > tmp && mv tmp asCHM.methyl_merged.txt
	cat title.txt asCHM.H3K9me3_merged.txt > tmp && mv tmp asCHM.H3K9me3_merged.txt

	# Score according to whether FC>=2 in both DNA methylation and H3K9me3
	# CHM with DNA methylation in maternal >= 0.5, H3K9me3 RPM in maternal >= 0.3,  fold change of DNA methylation level and H3K9me3 RPM (maternal / paternal)  >= 2 scored 1 point;
	# CHM with DNA methylation in paternal >= 0.5, H3K9me3 RPM in paternal >= 0.3, fold change of DNA methylation and H3K9me3 RPM (paternal / maternal) >= 2 scored -1 point;
	# CHM in other conditions scored 0 point.
	# Then we summed the DNA methylation and H3K9me3 asymmetry score in the five stages as its epigenetic modification asymmetry

	paste asCHM.methyl_merged.txt asCHM.H3K9me3_merged.txt | tail -n +2 | awk '
	 function score_E(M_mat, M_pat, H_mat, H_pat){RES=0;if(M_mat!="NA" && M_pat!="NA" && H_mat!="NA" && H_pat!="NA"){\
	 if(M_mat>=0.5 && H_mat>=0.3){if(M_pat/M_mat<=0.5 && H_pat/H_mat<=0.5){RES=1;}} \
	 else if(M_pat>=0.5 && H_pat>=0.3){if(M_mat/M_pat<=0.5 && H_mat/H_pat<=0.5){RES=-1;}}};return RES};\
	 BEGIN{FS=OFS="\t"}{print $1, $2, $3, $4, score_E($5, $6, $19, $20)+score_E($7, $8, $21, $22)+score_E($9, $10, $23, $24)+score_E($11, $12, $25, $26)+score_E($13, $14, $27, $28)}' > asCHM.score_EpiAsy.txt
}
epiAsym_F

expr_F(){
	# gene/transposable element with FPKM in GG >= 3 and fold change (GG/AG) >= 2 scored 1 point;
	# gene/transposable element with FPKM in AG >= 3 and fold change (AG/GG) >= 2 scored -1 point;
	# gene/transposable element in other conditions scored 0 point.
	# Then we summed the expression asymmetry score in the four stages, and was denoted as its summarized expression asymmetry score

	knownIGs_F(){
		windowBed -w 300000 -c -a asCHM.bed -b Expression/imprinted_genes_merged.promoter.bed | awk 'BEGIN{FS=OFS="\t"}{if($5>0){print $1, $2, $3, $4, 1} else{print $1, $2, $3, $4, 0}}' > asCHM.within300kb.knownIGs.txt
	}

	gene_F(){
		echo -e "#Chrom\tStart\tEnd\tName\tGenesWithin300kb\tAlleleScoreOfGenesWithin300kb\tAlleleScoreSummarized" > title.asCHM.score_expr.txt
		nTotal=$(cat asCHM.bed | wc -l)
		for nRow in $(seq 1 ${nTotal})
		do
			sed -n "${nRow}p" asCHM.bed > nEle.bed # temporary

			# Genes within 300kb
			windowBed -w 300000 -a nEle.bed -b ${dirPATH}/Expression/mm10.refGene.exprPattern.txt > nEle.AlleleScore.txt # temporary

			# Score of genes within 300kb
			nGenes=$(cat nEle.AlleleScore.txt | awk 'BEGIN{FS="\t";ORS=";"}{print $8"["$9"]"}')
			nGenesScore=$(cat nEle.AlleleScore.txt | awk 'function abs(v){return v > 0 ? v : -1*v};BEGIN{FS="\t";ORS=";"}{if(abs($11+$12+$13+$14)>=2){print $11+$12+$13+$14} else{print 0}}')
			nScore=$(cat nEle.AlleleScore.txt | awk 'function abs(v){return v > 0 ? v : -1*v};BEGIN{FS=OFS="\t";M=0;}{if(abs($11+$12+$13+$14)>abs(M) && abs($11+$12+$13+$14)>=2){M=$11+$12+$13+$14}}END{print M}')
			echo -e "${nGenes}\t${nGenesScore}\t${nScore}" > nEle.AlleleScoreManaged.txt # temporary

			paste nEle.bed nEle.AlleleScoreManaged.txt >> asCHM.score_geneAsy.txt
			rm nEle.bed nEle.AlleleScore.txt nEle.AlleleScoreManaged.txt
			unset nGenes nGenesScore nScore
		done # for nRow end
	}
	gene_F

	transposableElement_F(){
		nTotal=$(cat asCHM.bed | wc -l)
		for nRow in $(seq 1 ${nTotal})
		do
			sed -n "${nRow}p" asCHM.bed > nEle.bed # temporary

			# TEs within 300kb
			windowBed -w 300000 -a nEle.bed -b mm10.te_tx.exprPattern.txt > nEle.AlleleScore.txt # temporary

			# Score of TEs within 300kb
			nGenes=$(cat nEle.AlleleScore.txt | awk 'BEGIN{FS="\t";ORS=";"}{print $8"["$9"]"}')
			nGenesScore=$(cat nEle.AlleleScore.txt | awk 'function abs(v){return v > 0 ? v : -1*v};BEGIN{FS="\t";ORS=";"}{if(abs($11+$12+$13+$14)>=2){print $11+$12+$13+$14} else{print 0}}')
			nScore=$(cat nEle.AlleleScore.txt | awk 'function abs(v){return v > 0 ? v : -1*v};BEGIN{FS=OFS="\t";M=0;}{if(abs($11+$12+$13+$14)>abs(M) && abs($11+$12+$13+$14)>=2){M=$11+$12+$13+$14}}END{print M}')
			echo -e "${nGenes}\t${nGenesScore}\t${nScore}" > nEle.AlleleScoreManaged.txt # temporary

			paste nEle.bed nEle.AlleleScoreManaged.txt >> asCHM.score_TEAsy.txt
			rm nEle.bed nEle.AlleleScore.txt nEle.AlleleScoreManaged.txt
			unset nGenes nGenesScore nScore
		done # for nRow end
	}
	transposableElement_F
}
expr_F

regulatoryElement_F(){
	methylatedZFP57MotifAsy_F(){
		nTotal=$(cat asCHM.bed | wc -l)
		for nRow in $(seq 1 ${nTotal})
		do
			sed -n "${nRow}p" asCHM.bed > nEle.bed # temporary

			# Overlap with asyMethylated ZFP57 motif
			intersectBed -wao -a nEle.bed -b RegulatoryElement/methylatedZFP57MotifAsy.txt > nEle.AlleleScore.txt # temporary

			# Score
			# No motif
			noMotif=$(cut -f 5 nEle.AlleleScore.txt | sort -u)
			if [ ${noMotif} == "." ]
			then
				EachScore=0
				Score=0					
			else
				EachScore=$(cat nEle.AlleleScore.txt | awk 'BEGIN{FS="\t";ORS=";"}{print $18}')
				Score=$(cat nEle.AlleleScore.txt | awk 'function abs(v){return v > 0 ? v : -1*v};BEGIN{FS=OFS="\t";M=0;}{if(abs($18)>abs(M)){M=$18}}END{print M}' )
			fi
			echo -e "${EachScore}\t${Score}" > nEle.AlleleScoreManaged.txt # temporary
			paste nEle.bed nEle.AlleleScoreManaged.txt >> asCHM.score_methylatedZFP57Asy.txt
			rm nEle.bed nEle.AlleleScore.txt nEle.AlleleScoreManaged.txt
			unset EachScore Score
		done # for nRow end
	}
	methylatedZFP57MotifAsy_F

	CTCFbinding_F(){
		intersectBed -c -a asCHM.bed -b RegulatoryElement/CTCF.bed | awk 'BEGIN{FS=OFS="\t"}{if($5==0){print $1, $2, $3, $4, 0} else{print $1, $2, $3, $4, 1}}' > asCHM.score_CTCFBinding.txt # Data source: Cistrome db
	}
	CTCFbinding_F

	lncRNAOverlap_F(){
		intersectBed -c -a asCHM.bed -b mm10.lncRNA_transcript.bed | awk 'BEGIN{FS=OFS="\t"}{if($5==0){print $1, $2, $3, $4, 0} else{print $1, $2, $3, $4, 1}}' > asCHM.score_lncRNAOverlap.txt # Data source: GENCODE
	}
	lncRNAOverlap_F
}
regulatoryElement_F

score_F(){
	echo -ne "#Chrom\tStart\tEnd\tName" > asCHM.score.title.txt
	cat asCHM > asCHM.score.txt

	# Asymmetric epigenetic modifications
	echo -ne "\tasym.epigenetic\tasym.modification" >> asCHM.score.title.txt
	cat asCHM.score_EpiAsy.txt | awk 'function abs(v){return v > 0 ? v : -1*v};BEGIN{FS=OFS="\t";}{print $5, abs($5)/5}' | \
	 paste asCHM.score.txt - > tmp && mv tmp asCHM.score.txt

	# Allele-specific expressed genes or TEs
	echo -ne "\tknown.IG\tasym.gene\tasym.TE\tasym.expression" >> asCHM.score.title.txt
	paste asCHM.within300kb.knownIGs.txt asCHM.score_geneAsy.txt asCHM.score_TEAsy.txt | awk 'function abs(v){return v > 0 ? v : -1*v};function max(a, b){return a > b ? a : b};BEGIN{FS=OFS="\t"}{print $5, $12/4, $19/4, max($5, max(abs($12), abs($19))/4)}' | \
	 paste asCHM.score.txt - > tmp && mv tmp asCHM.score.txt # MD5 is same

	# Regulatory element
	echo -ne "\tasym.methylated.ZFP57\tCTCForlncRNA\tregulatory.elements" >> asCHM.score.title.txt
	paste asCHM.score_methylatedZFP57Asy.txt asCHM.score_CTCFBinding.txt asCHM.score_lncRNAOverlap.txt | awk 'function abs(v){return v > 0 ? v : -1*v};function max(a, b){return a > b ? a : b};BEGIN{FS=OFS="\t"}{print $6/5, max($11, $16), (abs($6/5)+max($11, $16))/2}' | \
	 paste asCHM.score.txt - > tmp && mv tmp asCHM.score.txt # MD5 is same 

	echo -ne "\tscore\n" >> asCHM.score.title.txt
	cat asCHM.score.txt | awk 'BEGIN{FS=OFS="\t"}{if($6>=0.2){print $0, $6+$10+$13} else{print $0, 0}}' | sort -k14,14nr -k7,7nr -k1,1 > tmp && mv tmp asCHM.score.txt

	cat asCHM.score.title.txt asCHM.score.txt > tmp && mv tmp asCHM.score.txt
}
score_F