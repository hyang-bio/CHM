import os, sys
import re
import twobitreader

USAGE = "USAGE: %prog <sequence> <bedFile> <outputFile> <2bit>"

def getSeqLoc(seq, bedF, outF, genomeF):
	bedFH = open(bedF, 'r')
	outFH = open(outF, 'w')
	genome = twobitreader.TwoBitFile(genomeF)
	for line in bedFH:
		ele = line.strip().split('\t')
		chrom = ele[0]
		start = ele[1]
		end = ele[2]

		sequences = genome[chrom][int(start):int(end)].upper()
		patterns = re.compile(seq)
		if seq in sequences:
			match = patterns.finditer(sequences)
			for s in match:
				loci = [int(start) + int(s.start()), int(start) + int(s.start()) + len(seq)]
				outFH.write("{0}\t{1}\t{2}\t{3}\t{4}\n" . format(chrom, loci[0], loci[1], start, end))
	outFH.close()
	bedFH.close()

def main():
	seq = sys.argv[1]
	bedF = sys.argv[2]
	outF = sys.argv[3]
	genomeF = sys.argv[4]
	getSeqLoc(seq, bedF, outF, genomeF)

main()