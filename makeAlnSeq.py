import sys, csv, subprocess, os, argparse

parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python 8-7-flair-to-fusions-pipe.py -i flair.aligned.bam -o outputPrefix -b buffer -a path to annotations')
parser.add_argument('-f', '--fastq', action='store', dest='f', default="", help='.fastq file')
parser.add_argument('-r', '--reads', action='store', dest='r', default="", help='.bed file')
args = parser.parse_args()

bedfile = open(args.r, 'r')
outfile = open(args.r.rstrip('.bed') + '-seq.bed', 'w')
genome = open(args.f, 'r')
stringlen = []
genomeDict = {}
lastKey = ""
for line in genome:
	if line[0][0] == '>':
		lastKey = line.lstrip('>').rstrip()
	else:
		genomeDict[lastKey] = line.rstrip().upper()

for line in bedfile:
	line = line.rstrip().split('\t')
	cigar, thisSeq = line[7], line[11]
	chrSize = len(genomeDict[line[0]])
	currString = "-" * int(line[1])
	if thisSeq != '*':
		while len(cigar) >= 2:
				i = 0
				while cigar[i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']:
					i += 1
				if cigar[i] in ['M', 'I', 'X', 'S']:
					if cigar[i] == 'M' or cigar[i] == 'X':
						currString += thisSeq[:int(cigar[:i])]
					thisSeq = thisSeq[int(cigar[:i]):]
				if cigar[i] == 'D' or cigar[i] == 'N':
					currString += '-' * int(cigar[:i])
				cigar = cigar[i+1:]
	else:
		thisChr = genomeDict[line[0]]
		thisSeq = thisChr[int(line[1]):]
		while len(cigar) >= 2:
			i = 0
			while cigar[i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']:
				i += 1
			if cigar[i] in ['M', 'X', 'D', 'N']:
				if cigar[i] == 'M':
					currString += thisSeq[:int(cigar[:i])]
				elif cigar[i] == 'X':
					currString += 'N' * int(cigar[:i])
				elif cigar[i] == 'D' or cigar[i] == 'N':
					currString += '-' * int(cigar[:i])
				thisSeq = thisSeq[int(cigar[:i]):]
			cigar = cigar[i+1:]
	currString += "-" * (chrSize-int(line[2]))
	line[11] = currString
	stringlen.append(len(currString))
	outfile.write('\t'.join(line) + '\n')
#print(sum(stringlen)/len(stringlen))
for chr in genomeDict:
	line = [chr.lstrip('>'), '0', str(len(genomeDict[chr])),'reference', '60', '+', '0', '*', '*', '0', '0', genomeDict[chr]]
	outfile.write('\t'.join(line) + '\n')
