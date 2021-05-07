# Main codes for  "Functional allele-specific H3K9me3 and DNA methylation co-marked CpG-rich regions in pre-implantation embryo".

##### 1. Identification of CpG-rich regomic loci with high H3K9me3 signal and DNA methylation level (CHM)

> To identify CpG-rich genomic loci with high H3K9me3 signal and DNA methylation (CHM) in each stage, we ran ChromHMM (v1.22) at the default 200-bp resolution in GG, AG and hybrid pre-implantation embryos, and 500-bp resolution in SNP-trackable E6.5 embryos as their relative low information. 

+ The source code for identification of CHM is available at [CHM_identification](CHM_identification).
---


##### 2. Identification of ICR-like regions (ICRLR) 
> To identify ICRLR with transcriptional regulation functions similar to those of known ICR, we calculated scores for all identified allele-specific CHMs based on features of known ICRs. Three groups of features were used, including allele-specific epigenetics, allele-specific expression and regulatory elements. 

+ The source code for identification of ICRLR is available at [ICRLR_identification](ICRLR_identification).
---