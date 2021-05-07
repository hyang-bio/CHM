#!/bin/bash

# ———————————————————————————————————
# To identify CpG-rich genomic loci with high H3K9me3 signal and DNA methylation (CHM) in each stage, 
# we ran ChromHMM (v1.22) at the default 200-bp resolution in gametes, GG, AG and hybrid pre-implantation embryos, 
# and 500-bp resolution in SNP-trackable E6.5 embryos as their relative low information.
# Take 200-bp resolution as an example here.

# mm10.2bit is downloaded from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit).
# mm10_euch.chrom.sizes is downloaded from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes) and only focus on those located in known euchromosome.
# mm10.b200bp.euchr.bed is annotations of 200-bp bins managed from mm10_euch.chrom.sizes.

# ${sample}.H3K9me3.rmDup.bam 
# ${sample}.sam.G.bed is output from mcall.
# ———————————————————————————————————


dirPATH=$(pwd)

# step1. Binarize H3K9me3 signal from bam file
cd ${dirPATH}
echo -e "${sample}\tH3K9me3\t${sample}.H3K9me3.rmDup.bam" > config_H3K9me3.txt
Java -Xmx50000m -jar ChromHMM.jar BinarizeBam -b 200 -paired mm10_euch.chrom.sizes ${dirPATH} config_H3K9me3.txt Binarized_H3K9me3

# step2. Binarize DNA methylation level from sam.G.bed in 200-bp bins (threshold = 0.5)
bash averageMethylInRegionMultipleThreads.sh mm10.b200bp.euchr.bed ${sample}.sam.G.bed ${sample}.methyl.txt 10
mkdir -p ${dirPATH}/Binarized_Methyl;cd ${dirPATH}/Binarized_Methyl
for chrom in $(seq 1 19)
do
	echo -e "${sample}\tchr${chrom}\nMethyl" > ${sample}_chr${chrom}_binary.txt
	grep -w "chr${chrom}" ../${sample}.signal_Methyl.txt | awk 'BEGIN{FS=OFS="\t"}{\
	 if($4=="NA"){print "2"} \
	 else if($4<0.5){print "0"} \
	 else{print "1"}}' >> ${sample}_chr${chrom}_binary.txt
done # for chrom end

# step3. Binarize CpG number in 200-bp bins (threshold = 6)
cd ${dirPATH}
python getCpGnumberMultipleTs.py mm10.b200bp.euchr.bed mm10.b200bp.euchr mm10.2bit
mkdir -p ${dirPATH}/Binarized_CpGNumber;cd ${dirPATH}/Binarized_CpGNumber
for chrom in $(seq 1 19)
do
	echo -e "200bp\tchr${chrom}\nCpGNumber" > 200bp_chr${chrom}_binary.txt
	for chrom in $(seq 1 19)
	do
		echo -e "200bp\tchr${chrom}\nCpGNumber" > 200bp_chr${chrom}_binary.txt
		grep -w "chr${chrom}" ${dirPATH}/mm10.b200bp.euchr.CpGnumber | awk 'BEGIN{FS=OFS="\t"}{\
		 if($4<6){print "0"} else{print "1"}}' >> 200bp_chr${chrom}_binary.txt
	done # for chrom end
done # for chrom end

# step4. Binarize according to co-occupancy of CpG number, H3K9me3 and DNA methylation
mkdir -p ${dirPATH}/Binarized_Cooccupancy;cd ${dirPATH}/Binarized_Cooccupancy
for chrom in $(seq 1 19)
do
	echo -e "${sample}\tchr${chrom}\nMerged" > ${sample}_chr${chrom}_binary.txt
	paste ${dirPATH}/Binarized_Methyl${methyl}/${sample}_chr${chrom}_binary.txt ${dirPATH}/Binarized_H3K9me3_TreatOnlyPValue${k9}/${sample}_chr${chrom}_binary.txt ${dirPATH}/Binarized_CpGNumber/200bp_chr${chrom}_binary.txt | \
	 tail -n +3 | awk 'BEGIN{FS=OFS="\t"}{if($1==2){print 2} else{if($1==1 && $2==1 && $3==1){print 1} else{print 0}}}' >> ${sample}_chr${chrom}_binary.txt
done # for chrom

# step5. Learn model
cd ${dirPATH}
java -Xmx50000m -jar ChromHMM.jar LearnModel -init random -p 20 Binarized_Cooccupancy Output 2 mm10
cd Output
State=$(cat emissions_2_999.txt | awk 'BEGIN{FS=OFS="\t";S="E1";MAX=-1;}{if($2+$3+$4>MAX){MAX=$2+$3+$4;S=$1}}END{print "E"S}')
grep -w ${State} ${sample}_2_999_segments.bed | cut -f 1-3 | sort -k1,1 -k2,2n | mergeBed -i - -d 2000 | awk 'BEGIN{FS=OFS="\t"}{if($3-$2>=600){print }}' > ${sample}.CHM.bed