#!/bin/bash

# ———————————————————————————————————
# Calculate average DNA methylation level
# ———————————————————————————————————

bedF=${1}
cpgF=${2}
outF=${3}
threads=${4}

# step1. Split into multiple small files
n_row=`cat ${bedF} | wc -l`
para_l=`bc <<< ${n_row}/${threads}+1`
split -l ${para_l} ${bedF} -d -a 3 ${outF}_subfile_

# step2. Calculate average DNA methylation for each small file
for file in `ls ${outF}_subfile_*`
do
	grep '^chr' ${file} | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3}' > ${file}.tmp && \
	intersectBed -wao -a ${file}.tmp -b ${cgF} > ${file}.all.G && \
	python methyl_processing.py region_methylation ${file}.all.G ${file}.ave &
done # for file end
wait

# step3. Merge into a single file
cat ${outF}_subfile_*.ave | sort -S100G -k1,1 -k2,2n --parallel=24 > ${outF}
rm ${outF}_subfile_*