import sys, csv, subprocess, os, argparse
from datetime import date
from statistics import median
from collections import Counter
import copy
import math
# import time
# import json
import numpy

def binarySearch(arr, t):
	if t <= arr[0]: return arr[0]
	if t >= arr[-1]: return arr[-1]
	i, j, mid = 0, len(arr)-1, 0
	while i < j:
		mid = int((i + j)/2)
		if arr[mid] == t: return arr[mid]
		elif t < arr[mid]:
			if mid > 0 and t > arr[mid-1]:
				if abs(arr[mid]-t) < abs(arr[mid-1]-t): return arr[mid]
				else: return arr[mid-1]
			j = mid
		else:
			if mid < len(arr)-1 and t < arr[mid+1]:
				if abs(arr[mid]-t) < abs(arr[mid+1]-t): return arr[mid]
				else: return arr[mid+1]
			i = mid + 1

def grouper(iterable, key):
    prev = None
    group = []
    for item in iterable:
        if prev is None or item[key] - prev[key] <= 15:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group

#Removing all repetitive region references

parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python 8-7-flair-to-fusions-pipe.py -i flair.aligned.bam -o outputPrefix -b buffer -a path to annotations')
parser.add_argument('-o', '--output', action='store', dest='o', default=date.today().strftime("%d-%m-%Y"), help='output file name base (default: date)')
parser.add_argument('-r', '--reads', action='store', dest='r', default="", help='.fa or fq file')
parser.add_argument('-m', '--bedFile', action='store', dest='m', default="", help='.bed file')
parser.add_argument('-f', '--flair', action='store', dest='f', default="/private/home/cafelton/flair-new/flair.py", help='flair path')
parser.add_argument('-g', '--genome', action='store', dest='g', default="/private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa", help='path to genome')
#parser.add_argument('-x', '--minimap', action='store', dest='x', default="/private/groups/brookslab/bin/minimap2", help='path to minimap')
parser.add_argument('-k', '--remapSize', action='store', dest='k', default=0, type=int, help='size of area to remap - only remaps if this is specified')
parser.add_argument('-t', '--transcriptome', action='store', dest='t', default="/private/groups/brookslab/reference_annotations/gencode.v37.annotation.gtf", help='path to transcriptome (.gtf)')
parser.add_argument('-n', '--spliceJunctions', action='store', dest='n', default="/private/groups/brookslab/cafelton/fusions-code/intropolis.liftover.hg38.junctions.sorted.txt", help='path to splice junction file (.txt)')
parser.add_argument('-b', '--buffer', action='store', dest='b', default=50000, help='length of buffer for combining nearby regions')
parser.add_argument('-l', '--readSupport', action='store', dest='l', default=3, help='number of reads required to call fusion')
parser.add_argument('-a', '--anno', action='store', dest='a', default=os.path.dirname(__file__) + '/gencode.v37.annotation-short.gtf', help='path to anno.gtf')
parser.add_argument('-p', '--bedProcess', action='store_true', dest='p', help='whether to take .bam and convert to .bed and process (True = assume existing processed .bam)')
parser.add_argument('-s', '--samConvert', action='store_true', dest='s', help='whether to convert .bam to .sam or (True = convert .bam (from fq prefix) to .sam)')
parser.add_argument('-y', '--includeMito', action='store_true', dest='y', help='whether to include fusions that are in the mitochondria (True=include)')
# parser.add_argument('-q', '--geneCov', action='store_true', dest='q', help='whether to filter out fusions ')
# parser.add_argument('-v', '--fastqCov', action='store_true', dest='v', help='whether to include fusions that are in the mitochondria (True=include)')
parser.add_argument('-w', '--removePromiscuous', action='store_true', dest='w', help='whether to filter out promiscuous chimeric genes')
parser.add_argument('-u', '--flairAlign', action='store_true', dest='u', help='whether to run flair align (True=already aligned, dont run)')
parser.add_argument('-c', '--flairCorrect', action='store_true', dest='c', help='whether to run flair correct (True=already corrected, dont run)')
parser.add_argument('-d', '--detectFusions', action='store_true', dest='d', help='whether to detect fusions (True=already detected, dont run)')
parser.add_argument('-i', '--callIsoforms', action='store_true', dest='i', help='whether to detect fusion isoforms (True=already detected, dont run)')
parser.add_argument('-j', '--matchFusionIsos', action='store_true', dest='j', help='whether to match isoforms to fusions (True=already matched or dont want to match, dont run)')
#/private/groups/brookslab/reference_annotations/
args = parser.parse_args()
# print('buffer', args.b)
prefix = '.'.join(args.r.split('.')[:-1])
bedFileMade = False
if len(args.r) > 0 and not args.u and not args.p:
	process = subprocess.Popen('python3 ' + args.f + ' align -g ' + args.g + ' -r ' + args.r + ' -o ' + prefix + '.aligned' +  #+ ' -m ' + args.x,
								'; bamToBed -bed12 -i ' + prefix + '.aligned.bam > ' + prefix + '.aligned.bed'
								'; rm ' + prefix + '.aligned.bam ' + prefix + '.aligned.bam.bai',
							   stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())
	bedFileMade = True

if (bedFileMade or os.path.exists(prefix + '.aligned.bed') or len(args.m) > 0) and not args.c and not args.p:
	correctQ = args.m if len(args.m) > 0 else prefix + '.aligned.bed'
	process = subprocess.Popen('python3 ' + args.f + ' correct -g ' + args.g + ' -f ' + args.t + ' -q ' + correctQ + ' -o ' + prefix,
							   stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())
if not args.p:
	correctQ = args.m if args.c else prefix + '_all_corrected.bed'
	#'python ' + os.path.dirname(os.path.realpath(__file__)) + '/standardizeBed.py' + ' -i ' + correctQ,
	myCommands = ['bedtools intersect -wao -a ' + correctQ.rstrip('bed').rstrip('.') + '.bed' + ' -b ' + args.a + ' > ' + prefix + '-bedtools-genes.txt',
				  'python ' + os.path.dirname(os.path.realpath(__file__)) + '/bedtoolsGeneHelper.py' + ' -i ' + prefix + '-bedtools-genes.txt',
				  'rm ' + prefix + '-bedtools-genes.txt']
	process = subprocess.Popen('; '.join(myCommands), stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())
	print('done with preprocessing')
filepath = '/'.join(prefix.split('/')[:-1]) + '/' if '/' in prefix else ''
outfilename = filepath + args.o + prefix.split('/')[-1]
if not args.d:
	args.b = int(args.b)
	meta = open(outfilename + "-meta.txt", "w")
	metadata = []
	bed = open(prefix + '-bedtools-genes-short.bed', 'r')
	junctions = {}
	count = 0
	readLength = {}
	last = None
	c = 0
	print('loading splice junctions and genes and paralogs')
	for line in open(args.n, 'r'):  # "/private/groups/brookslab/cafelton/fusions-code/gencode.v37.junctions.txt", 'r'):
		count += 1
		if count % 2 > 0:
			last = line.strip()
		else:
			temp = numpy.fromstring(line.rstrip(), dtype=numpy.int, sep=",")
			junctions[last] = temp
	allGeneLoc = {}
	for line in open(args.t, 'r'):
		if line[0] != "#":
			line = line.split('\t')
			if line[2] == 'gene':
				name = line[8].split('gene_name "')[1].split('"')[0]
				allGeneLoc[name] = [[line[0], int(line[3]), int(line[4]), line[6]]]
			elif line[2] == 'transcript':
				name = line[8].split('gene_name "')[1].split('"')[0]
				allGeneLoc[name].append([line[0], int(line[3]), int(line[4]), line[6]])
	paralogs = {}
	for line in open(os.path.dirname(os.path.realpath(__file__)) + '/paralog_clusters_with_dup.txt'):
		line = line.rstrip().split('\t')
		for i in line: paralogs[i] = line
	clinicalF = []
	for line in open("/private/groups/brookslab/cafelton/fusions-code/treehouse-clinical-fusions.txt"):
		clinicalF.append(frozenset(line.strip().split('--')))
	clinicalF = set(clinicalF)
	potential_chimeric = {}  # {read name: [entries]}
	print("finding chimeric reads")
	bedLines = []
	buffer = args.b
	maxMapQ, bedLineCount, avgMapQ = 0, 0, 0
	c = 0
	geneReads = {}
	#GET POTEINTIAL CHIMERIC (MULTIPLE MAPPING) READS FROM BED FILE
	c = 0
	for line in bed:
		if len(line) > 20:
			c += 1
			bedLineCount += 1
			line = line.rstrip().split('\t')
			bedLines.append(line)
			readname, gene  = line[3].split('--')
			#geneReads is important for understanding what fraction of reads at a locus are in the fusion
			if gene.split('/')[0] not in geneReads.keys():
				geneReads[gene.split('/')[0]] = []
			geneReads[gene.split('/')[0]].append(readname)
			if '/' not in gene: gene = gene.replace('.0', '')
			else: gene = gene.split('/')[0]
			avgMapQ += int(line[4])
			if int(line[4]) > maxMapQ: maxMapQ = int(line[4])
			if readname in potential_chimeric:
				if gene in potential_chimeric[readname]: potential_chimeric[readname][gene].append(line)
				else: potential_chimeric[readname][gene] = [line]
			else: potential_chimeric[readname] = {gene: [line]}

	fusions_found = {}  # {fused genes: count}
	fusionReads = []
	avgMapQ = avgMapQ/float(bedLineCount)
	print("finding potential fusions")
	c = 0
	for read in potential_chimeric:
		if len(potential_chimeric[read]) > 1:
			c += 1
			locs = list(potential_chimeric[read].keys())
			locs.sort()
			fusion_name = frozenset(locs)
			if fusion_name not in fusions_found:
				fusions_found[fusion_name] = {'mapScore':0, 'readNames':[]}
				for loc in locs:
					fusions_found[fusion_name][loc] = {'reads':[], 'left':[], 'right':[], 'strand':[], 'chr':[potential_chimeric[read][loc][0][0]]}
			for loc in locs:
				fusions_found[fusion_name][loc]['reads'] += potential_chimeric[read][loc] #This gives the whole line from the bed file
				for i in potential_chimeric[read][loc]:
					fusions_found[fusion_name]['mapScore'] += int(i[4])
					fusions_found[fusion_name][loc]['left'].append(int(i[1]))
					fusions_found[fusion_name][loc]['right'].append(int(i[2]))
					fusions_found[fusion_name][loc]['strand'].append(i[5])
			fusions_found[fusion_name]['readNames'].append(read)
	# print('fusions',  len(fusions_found.keys()), c)
	#AGGREGATE AND SORT NON-GENIC REGIONS IN FUSIONS
	print('condensing fusions in non-genic regions')
	locsToCondense = {}
	c = 0
	for i in fusions_found:
		locs = list(i)
		for j in locs:
			if len(j.split('-')) > 1 and j[:3] == 'chr':
				chr, loc = j.split('-')
				if chr not in locsToCondense: locsToCondense[chr] = []
				locsToCondense[chr].append(int(loc))
	#JOIN CLOSE REGIONS TOGETHER AND MARK THEM FOR UPDATING
	updateValues = {}
	for chr in locsToCondense:
		lastLoc = 0
		lastKey = None
		locsToCondense[chr].sort()
		for j in locsToCondense[chr]:
			if j-lastLoc <= buffer or j-lastLoc <= 1000 and lastKey == None:
				lastKey = '-'.join([chr, str(lastLoc)])
				updateValues[lastKey] = lastKey
				updateValues['-'.join([chr, str(j)])] = lastKey
			elif j-lastLoc <= buffer or j-lastLoc <= 1000:
				updateValues['-'.join([chr, str(j)])] = lastKey
			else:
				updateValues['-'.join([chr, str(j)])] = '-'.join([chr, str(j)])
				lastKey = None
			lastLoc = j

	#CONDENSE NON-GENIC FUSION REGIONS INTO FEWER FUSIONS
	new_fusions_found = {}
	# allMatches = []
	readToFusion = {}
	for i in fusions_found:
		locs, locsstatic = list(i), list(i)
		paralocs = {}
		for l in range(len(locs)):
			if len(locs[l].split('-')) > 1 and locs[l][:3] == 'chr':
				locs[l] = updateValues[locs[l]]
			if locs[l] in paralogs:
				if paralogs[locs[l]][0] not in paralocs: paralocs[paralogs[locs[l]][0]] = []
				paralocs[paralogs[locs[l]][0]].append(l)
			else: paralocs[locs[l]] = [l]
		newparalocs = {}
		for p in paralocs:
			a = []
			for l in paralocs[p]: a.append(locs[l])
			a.sort()
			if p not in a: newparalocs[a[0]] = paralocs[p]
			else: newparalocs[p] = paralocs[p]
		newlocs = list(newparalocs.keys())
		newlocs.sort()
		if len(newparalocs.keys()) <= 1:
			metadata.append(['--'.join(locsstatic), 'paralogs', 'pl', fusions_found[i]['readNames']])
		else:
			fusionname = frozenset(newlocs)
			if fusionname not in new_fusions_found.keys():
				new_fusions_found[fusionname] = {}
				new_fusions_found[fusionname]['mapScore'] = fusions_found[i]['mapScore']
				new_fusions_found[fusionname]['readNames'] = fusions_found[i]['readNames']
				# allMatches.extend(fusions_found[i]['readNames'])
				for r in fusions_found[i]['readNames']: readToFusion[r] = fusionname
				for l in newparalocs:
					new_fusions_found[fusionname][l] = fusions_found[i][locsstatic[newparalocs[l][0]]]
					for subl in newparalocs[l][1:]:
						for key in ['reads', 'left', 'right', 'strand', 'chr']:
							new_fusions_found[fusionname][l][key] += fusions_found[i][locsstatic[subl]][key]
			else:
				new_fusions_found[fusionname]['mapScore'] += fusions_found[i]['mapScore']
				new_fusions_found[fusionname]['readNames'].extend(fusions_found[i]['readNames'])
				# allMatches.extend(fusions_found[i]['readNames'])
				for r in fusions_found[i]['readNames']: readToFusion[r] = fusionname
				for l in newparalocs:
					if l not in new_fusions_found[fusionname]:
						new_fusions_found[fusionname][l] = fusions_found[i][locsstatic[newparalocs[l][0]]]
						for subl in newparalocs[l][1:]:
							for key in ['reads', 'left', 'right', 'strand', 'chr']:
								new_fusions_found[fusionname][l][key] += fusions_found[i][locsstatic[subl]][key]
					else:
						for subl in newparalocs[l]:
							for key in ['reads', 'left', 'right', 'strand', 'chr']:
								new_fusions_found[fusionname][l][key] += fusions_found[i][locsstatic[subl]][key]
	print('new fusions',  len(new_fusions_found.keys()))

	# print(new_fusions_found[frozenset(['SLC4A8','CAPN7'])])
	correctQ = prefix if os.path.exists(prefix + '.bam') else prefix + '.aligned'
	if args.s:
		print('making sam file')
		process = subprocess.Popen('samtools view -h -o ' + correctQ + '.sam ' + correctQ + '.bam',stdout=subprocess.PIPE, shell=True)
		print(process.communicate()[0].strip())
	correctQ = prefix if os.path.exists(prefix + '.sam') else prefix + '.aligned'
	# allMatches = set(allMatches)
	readLength = {}
	print('processing sam file')
	sam = open(correctQ + '.sam', 'r')
	flagbinary = {'0':0, '2048':0, '16':1, '2064':1}
	mappingLocs = {}
	fgenes = []
	finalchimeras = []
	for line in sam:
		line = line.rstrip().split('\t')
		if line[0] in readToFusion.keys():
			locs = []
			if readToFusion[line[0]] not in mappingLocs: mappingLocs[readToFusion[line[0]]] = {}
			if line[0] not in mappingLocs[readToFusion[line[0]]]: mappingLocs[readToFusion[line[0]]][line[0]] = []
			i = 0
			while line[5][i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']: i += 1
			if line[5][i] == 'M': locs.append(0)
			else: locs.append(int(line[5][:i]))
			i = -2
			while line[5][i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']: i -= 1
			if line[0] not in readLength: readLength[line[0]] = len(line[9])
			if line[5][-1] == 'M': locs.append(readLength[line[0]])
			else: locs.append(readLength[line[0]] - int(line[5][i + 1:-1]))
			currDist = 1000000000000000000000
			currGene = None
			for gene in list(readToFusion[line[0]]):
				loc = median(new_fusions_found[readToFusion[line[0]]][gene]['left'])
				locchr = max(set(new_fusions_found[readToFusion[line[0]]][gene]['chr']),
								 key=new_fusions_found[readToFusion[line[0]]][gene]['chr'].count)
				if locchr == line[2] or currGene == None:
					if abs(int(line[3]) - loc) < currDist:
						currDist, currGene = abs(int(line[3]) - loc), gene
			index = new_fusions_found[readToFusion[line[0]]]['readNames'].index(line[0])
			if flagbinary[line[1]] == 1:
				temp = locs[0]
				locs[0] = readLength[line[0]] - locs[1]
				locs[1] = readLength[line[0]] - temp
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[0], new_fusions_found[readToFusion[line[0]]][currGene]['right'][index], currGene])
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[1], new_fusions_found[readToFusion[line[0]]][currGene]['left'][index], currGene])
			else:
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[0], new_fusions_found[readToFusion[line[0]]][currGene]['left'][index], currGene])
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[1], new_fusions_found[readToFusion[line[0]]][currGene]['right'][index], currGene])




	print('filtering fusions', len(mappingLocs.keys()))
	filteredFusions, fusionlocbounds = [], {}
	longestfastqdist, readfastqdist = {}, {}
	readcorrectedbp = {}
	trueBreakpoints = {'ALK':[2, 29223529], 'EML4':[2,42264732], 'GUCY1A2':[1,106776438], 'PIWIL4':[1,94568729], 'SCAMP2':[1,74845472], 'WDR72':[2,53706075],
					   'RET':[1,43116583], 'CCDC6':[1,59906121], 'CHRNA6':[2,42765205], 'ERGIC2':[1,29366876], 'UBR3':[1,170061078], 'EFHD1':[2,232672444]}
	# tempout = open("EML4-breakpoint-extra-exon.tsv", 'w')
	# tempout.write('\t'.join(['dist of genomic mapping', 'mapping gap on read', 'genomic breakpoint location', 'extra exon?']) + '\n')
	# tempout = open("drr-fusion-loc-bounds.bed", 'w')
	tempout2 = open("drr-breakpoint-clusters.tsv", 'w')
	for f in mappingLocs:
		if len(mappingLocs[f].keys()) >= args.l:
			fastqcov, genecov, tooshort, genecovfinal, shortmap, isMito, mitogene = [], {}, False, [0], 0, False, None
			for l in list(f):
				new_fusions_found[f][l]['chr'] = max(set(new_fusions_found[f][l]['chr']), key=new_fusions_found[f][l]['chr'].count)
				genecov[l] = []
				for i in range(len(new_fusions_found[f][l]['left'])):
					genecov[l].append([new_fusions_found[f][l]['right'][i]-new_fusions_found[f][l]['left'][i], new_fusions_found[f][l]['left'][i], new_fusions_found[f][l]['right'][i]])
				if (not(l[:3] == 'chr' and '-' in l)) and l in allGeneLoc:
					if allGeneLoc[l][0][0] == 'chrM': isMito, mitogene = True, l
				elif l[:4] == 'chrM': isMito, mitogene = True, l
			if isMito:
				metadata.append(['--'.join(list(f)), 'isMito', mitogene, list(mappingLocs[f].keys())])
				continue
			for l in genecov:
				genecov[l].sort()
				# if median([x[0] for x in genecov[l][:int(len(genecov[l])/4)]]) < 200: tooshort = True
				midcov, bestcov = genecov[l][int(len(genecov[l])/2)], []
				if midcov[0] < 100: tooshort, shortmap = True, midcov[0]
				if not (l[:3] == 'chr' and '-' in l) and l in allGeneLoc:
					for t in allGeneLoc[l][1:]:
						if midcov[1] >= t[1]-50 and midcov[2] <= t[2]+50: bestcov.append(float(midcov[0])/(t[2]-t[1]))
					if len(bestcov) > 0: genecovfinal.append(sum(bestcov)/len(bestcov))
					else: genecovfinal.append(0)
			if max(genecovfinal) > 0.95:
				metadata.append(['--'.join(list(f)), 'geneCov', max(genecovfinal), list(mappingLocs[f].keys())])
				continue
			if tooshort:
				metadata.append(['--'.join(list(f)), 'tooShortMapping', shortmap, list(mappingLocs[f].keys())])
				continue
			minlen, alllen = 100000000, []
			for r in mappingLocs[f]:
				if len(mappingLocs[f][r]) < minlen: minlen = len(mappingLocs[f][r])
			possPromLocs = [[] for i in range(minlen)]
			firstgene = max(set([mappingLocs[f][x][0][2] for x in mappingLocs[f]]), key=[mappingLocs[f][x][0][2] for x in mappingLocs[f]].count)
			breakpoints = {}
			ends = {}
			for r in mappingLocs[f]:
				#mappingLocs[fusion][read] = [[loc on fastq, loc on genome, gene], [etc]]
				tempfastqcov = 0
				for l in range(0, len(mappingLocs[f][r]), 2):
					tempfastqcov += abs(mappingLocs[f][r][l+1][0]-mappingLocs[f][r][l][0])
				fastqcov.append(float(tempfastqcov)/readLength[r])
				mappingLocs[f][r].sort()
				x = 1
				while x < len(mappingLocs[f][r])-2:
					#condensing adjacent mappings in the same gene
					if mappingLocs[f][r][x][2] == mappingLocs[f][r][x + 1][2]:
						mappingLocs[f][r] = mappingLocs[f][r][:x] + mappingLocs[f][r][x+2:]
					if mappingLocs[f][r][x][2] != mappingLocs[f][r][x-1][2]:
						temp = mappingLocs[f][r][x]
						mappingLocs[f][r][x] = mappingLocs[f][r][x+1]
						mappingLocs[f][r][x+1] = temp
					x += 2
				while len(mappingLocs[f][r]) > minlen:
					#removing extra mappings on the edges - if A-B-A, convert to A-B - seems like a bad strategy, should remove later
					if mappingLocs[f][r][0][2] == mappingLocs[f][r][-1][2]:
						if abs(mappingLocs[f][r][0][0]-mappingLocs[f][r][1][0]) < abs(mappingLocs[f][r][-1][0]-mappingLocs[f][r][-2][0]):
							mappingLocs[f][r] = mappingLocs[f][r][2:]
						else: mappingLocs[f][r] = mappingLocs[f][r][:-2]
					else: mappingLocs[f][r] = mappingLocs[f][r][:-2]
				#just making sure all reads are oriented in the same way - read 1: A-B, read 2: B-A -> read1:A-B, read2: A-B
				if mappingLocs[f][r][0][2] != firstgene:
					mappingLocs[f][r] = mappingLocs[f][r][::-1]
					for x in range(len(mappingLocs[f][r])): mappingLocs[f][r][x][0] = readLength[r] - mappingLocs[f][r][x][0]
				for x in range(1, len(mappingLocs[f][r])-2, 2):
					bpname, gap = frozenset([mappingLocs[f][r][x][2], mappingLocs[f][r][x + 1][2]]), abs(mappingLocs[f][r][x][0] - mappingLocs[f][r][x + 1][0])
					#this is just for my testing of fastq dist corr to mapping did to breakpoint, not essential
					if r in readfastqdist:
						if gap > readfastqdist[r]: readfastqdist[r] = gap
					else: readfastqdist[r] = gap
					if bpname not in breakpoints:
						breakpoints[bpname] = {'fastqdist': [], 'genomedist': []}
						for y in range(2): breakpoints[bpname][mappingLocs[f][r][x+y][2]] = []
					breakpoints[bpname]['fastqdist'].append(gap)
					if new_fusions_found[f][mappingLocs[f][r][x][2]]['chr'] == new_fusions_found[f][mappingLocs[f][r][x + 1][2]]['chr']:
						breakpoints[bpname]['genomedist'].append(abs(mappingLocs[f][r][x][1] - mappingLocs[f][r][x + 1][1]))
					breakpoints[bpname][mappingLocs[f][r][x][2]].append([abs(mappingLocs[f][r][x][1]-mappingLocs[f][r][x-1][1]), gap, mappingLocs[f][r][x][1]])
					breakpoints[bpname][mappingLocs[f][r][x+1][2]].append([abs(mappingLocs[f][r][x+1][1] - mappingLocs[f][r][x+2][1]), gap, mappingLocs[f][r][x+1][1]])
					if r not in readcorrectedbp: readcorrectedbp[r] = {}
					readcorrectedbp[r][mappingLocs[f][r][x][2]] = [mappingLocs[f][r][x][1], binarySearch(junctions[new_fusions_found[f][mappingLocs[f][r][x][2]]['chr']],mappingLocs[f][r][x][1])]
					readcorrectedbp[r][mappingLocs[f][r][x+1][2]] = [mappingLocs[f][r][x+1][1], binarySearch(junctions[new_fusions_found[f][mappingLocs[f][r][x+1][2]]['chr']], mappingLocs[f][r][x+1][1])]
				for x in [0, -1]:
					if mappingLocs[f][r][x][2] not in ends: ends[mappingLocs[f][r][x][2]] = []
				#dist of mapping on genome, dist of mapping on fastq, loc on genome
				ends[mappingLocs[f][r][0][2]].append([abs(mappingLocs[f][r][0][1]-mappingLocs[f][r][1][1]), abs(mappingLocs[f][r][0][0]-mappingLocs[f][r][1][0]), mappingLocs[f][r][0][1]])
				ends[mappingLocs[f][r][-1][2]].append([abs(mappingLocs[f][r][-1][1] - mappingLocs[f][r][-2][1]),abs(mappingLocs[f][r][-1][0] - mappingLocs[f][r][-2][0]),mappingLocs[f][r][-1][1]])
				for l in range(len(mappingLocs[f][r])): possPromLocs[l].append(mappingLocs[f][r][l])
			if median(fastqcov) < 0.8: metadata.append(['--'.join(list(f)), 'fastqCov', median(fastqcov), list(mappingLocs[f].keys())])
			# print("breakpoints", breakpoints)
			tooclose, fastqdistfinal, closedist, ssdist = False, 0, 0, []
			geneorder = []
			for bp in breakpoints:
				# print(f, breakpoints[bp]['genomedist'])
				if geneorder == []: geneorder = list(bp)
				else:
					temp = list(bp)
					print(geneorder, temp)
					for i in range(2):
						if temp[i] in geneorder:
							if geneorder.index(temp[i]) == 0: geneorder = [temp[abs(i-1)]] + geneorder
							else: geneorder.append(temp[abs(i-1)])
							break
				# print(bp, sorted(breakpoints[bp]['fastqdist']))
				if median(breakpoints[bp]['fastqdist']) > 15: fastqdistfinal = median(breakpoints[bp]['fastqdist'])
				if len(breakpoints[bp]['genomedist']) > 0 and median(breakpoints[bp]['genomedist']) < args.b: tooclose, closedist = True, median(breakpoints[bp]['genomedist'])
				for loc in set(breakpoints[bp].keys()) - {'fastqdist', 'genomedist'}:
					#just for coloring reads by fastq dist
					if loc in longestfastqdist:
						if max(breakpoints[bp]['fastqdist']) > longestfastqdist[loc]: longestfastqdist[loc] = max(breakpoints[bp]['fastqdist'])
					else: longestfastqdist[loc] = max(breakpoints[bp]['fastqdist'])

					breakpoints[bp][loc].sort(key=lambda x: x[2])
					c = 0
					for group in grouper(breakpoints[bp][loc], 2):
						if len(group) >= 0.1 * len(breakpoints[bp][loc]) and median([x[1] for x in group]) < 50:
							tempout2.write('\t'.join(['--'.join(list(f)), loc, str(c), str(len(group)), str(median([x[2] for x in group])), str(median([x[1] for x in group]))]) +'\n') #+ '\t' + ','.join([str(x[1]) + '-' +  str(x[2]) for x in group]) + '\n')
							c += 1
					breakpoints[bp][loc].sort()
					#selecting for longest mapping length on genome
					breakpoints[bp][loc] = breakpoints[bp][loc][-1*int(math.ceil(float(len(breakpoints[bp][loc]))/4)):][::-1]
					breakpoints[bp][loc].sort(key=lambda x:abs(x[1])) #sorting by shortest mapping gap on fastq
					# if loc == 'EML4':
					# 	for l in breakpoints[bp][loc]:
					# 		extraexon = 'extra' if abs(l[2]-42264730) > abs(l[2]-42264950) else 'correct'
					# 		tempout.write('\t'.join([str(x) for x in l] + [extraexon]) + '\n')
					closestSS = binarySearch(junctions[new_fusions_found[f][loc]['chr']],int(median([x[2] for x in breakpoints[bp][loc][:int(math.ceil(float(len(breakpoints[bp][loc]))/4))]])))
					# print(loc, int(median([x[2] for x in breakpoints[bp][loc][:int(math.ceil(float(len(breakpoints[bp][loc]))/4))]])), closestSS, trueBreakpoints[loc] - closestSS)
					ssdist.append(abs(closestSS - breakpoints[bp][loc][0][2]))
					breakpoints[bp][loc] = closestSS

					# if f == frozenset(['SLC25A24', 'NBPF6']): print(breakpoints[bp][loc])
			if fastqdistfinal > 15: metadata.append(['--'.join(list(f)), 'fastqDist', fastqdistfinal, list(mappingLocs[f].keys())])
			if tooclose: metadata.append(['--'.join(list(f)), 'tooClose', closedist, list(mappingLocs[f].keys())])
			if max(ssdist) > 10: metadata.append(['--'.join(list(f)), 'ssDist', max(ssdist), list(mappingLocs[f].keys())])
			distToProm, gene5 = 1000000000000000000000, None
			# if 'RET' in f: print("ends", ends)
			fusionlocbounds[f] = {}
			for loc in ends:
				ends[loc].sort(key=lambda x:abs(x[1]), reverse=True)
				ends[loc] = ends[loc][:int(math.ceil(float(len(ends[loc]))/4))]
				ends[loc].sort(reverse=True)
				medianEnd = int(median([x[2] for x in ends[loc][:int(math.ceil(float(len(ends[loc]))/4))]]))
				fusionlocbounds[f][loc] = [medianEnd]
				if loc in allGeneLoc:
					for t in allGeneLoc[loc][1:]:
						genestart = t[1] if t[3] == '+' else t[2]
						if abs(medianEnd - genestart) < distToProm:
							distToProm = abs(medianEnd - genestart)
							gene5 = loc
			# print(f, median(fastqcov), median(genecovfinal), fastqdistfinal, tooclose, ssdist, distToProm, gene5) #if f == frozenset(['SLC25A24', 'NBPF6']):
			if distToProm > 300: metadata.append(['--'.join(list(f)), 'distToProm', distToProm, list(mappingLocs[f].keys())])
			if gene5 != geneorder[0]: geneorder = geneorder[::-1]
			#USE THIS FOR ALL FILES EXCEPT DRR and distToProm <= 300 and fastqdistfinal <= 15
			if not isMito and median(fastqcov) >= 0.8 and max(ssdist) <= 10 and median(genecovfinal) <= 0.95 and not tooshort and not tooclose:
				for g in f:
					fgenes.append(g)
				finalchimeras.append(['--'.join(list(f)), list(mappingLocs[f].keys())])
				locbreakpoints = []
				for i in range(0, len(geneorder)-1):
					locbreakpoints.append([geneorder[i], new_fusions_found[f][geneorder[i]]['chr'], str(breakpoints[frozenset([geneorder[i], geneorder[i+1]])][geneorder[i]])])
					if geneorder[i] not in fusionlocbounds[f]: fusionlocbounds[f][geneorder[i]] = []
					fusionlocbounds[f][geneorder[i]].append(breakpoints[frozenset([geneorder[i], geneorder[i+1]])][geneorder[i]])
				locbreakpoints.append([geneorder[-1], new_fusions_found[f][geneorder[-1]]['chr'], str(breakpoints[frozenset([geneorder[-2], geneorder[-1]])][geneorder[-1]])])
				fusionlocbounds[f][geneorder[-1]].append(breakpoints[frozenset([geneorder[-1], geneorder[-2]])][geneorder[-1]])
				filteredFusions.append(['--'.join(geneorder), len(mappingLocs[f]), round(new_fusions_found[f]['mapScore']/float(len(mappingLocs[f])), 2)] + ['-'.join(x) for x in locbreakpoints])
				for loc in fusionlocbounds[f]:
					fusionlocbounds[f][loc].sort()
					# tempout.write('\t'.join([new_fusions_found[f][loc]['chr'], str(fusionlocbounds[f][loc][0]), str(fusionlocbounds[f][loc][1]), loc]) + '\n')
		else: metadata.append(['--'.join(list(f)), 'readSupport', len(mappingLocs[f].keys()), list(mappingLocs[f].keys())])
	for line in finalchimeras:
		thesegenes = line[0].split('--')
		if fgenes.count(thesegenes[0]) > 2 or fgenes.count(thesegenes[1]) > 2:
			metadata.append([line[0], 'promiscuous', str(fgenes.count(thesegenes[0])) + ',' + str(fgenes.count(thesegenes[1])), line[1]])
		else:
			metadata.append([line[0], 'predictedSV', str(fgenes.count(thesegenes[0])) + ',' + str(fgenes.count(thesegenes[1])),line[1]])
	reads = open(outfilename + "Reads.bed", "w")
	#reads.write('track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"\n')
	# print(fusionlocbounds)
	for i in metadata:
		i[3] = ','.join(i[3])
		meta.write('\t'.join([str(x) for x in i]) + '\n')#i + "\t" + str(metadata[i]) + "\n")
	meta.close()
	# print(longestfastqdist)
	# c = 0
	# for x in readcorrectedbp:
	# 	if c < 5: print(x, readcorrectedbp[x])
	# 	c += 1
	fastqdistvsbp = [[],[]]
	fusions = open(outfilename + "Fusions.tsv", "w")
	filteredFusions.sort(key=lambda x:x[1], reverse=True)
	c = 0
	for line in filteredFusions:
		fusions.write('\t'.join([str(x) for x in line]) + '\n')
		for l in list(set(line[0].split('--'))):
			# print(line[0].split('--'))
			# print([r[3] for r in new_fusions_found[frozenset(line[0].split('--'))][l]['reads']])
			for r in new_fusions_found[frozenset(line[0].split('--'))][l]['reads']:
				if r[3].split('--')[1].split('/')[0] in trueBreakpoints:
					# print(line)
					# print(r[3].split('--')[1].split('/')[0], trueBreakpoints[r[3].split('--')[1].split('/')[0]])
					# print(abs(int(r[trueBreakpoints[r[3].split('--')[1].split('/')[0]][0]])-trueBreakpoints[r[3].split('--')[1].split('/')[0]][1]), int(line[trueBreakpoints[r[3].split('--')[1].split('/')[0]][0]]))
					# if abs(int(r[trueBreakpoints[r[3].split('--')[1].split('/')[0]][0]])-trueBreakpoints[r[3].split('--')[1].split('/')[0]][1]) > 160000: print(line[0], r)
					fastqdistvsbp[0].append(readfastqdist[r[3].split('--')[0]])
					fastqdistvsbp[1].append(abs(int(r[trueBreakpoints[r[3].split('--')[1].split('/')[0]][0]])-trueBreakpoints[r[3].split('--')[1].split('/')[0]][1]))
				currGene = r[3].split('--')[1].split('/')[0]
				if int(r[1]) > fusionlocbounds[frozenset(line[0].split('--'))][currGene][0] - 150 and int(r[2]) < fusionlocbounds[frozenset(line[0].split('--'))][currGene][1] + 150:
					r[8] = str(int((250/math.log(longestfastqdist[r[3].split('--')[1].split('/')[0]] + 1))*math.log(readfastqdist[r[3].split('--')[0]] + 1) + 5)) + ",0,255"
					r[3] = '-.-'.join([line[0]] + r[3].split('--') + [str(readfastqdist[r[3].split('--')[0]])])
					# c += 1
					for x in [1,2]:
						if int(r[x]) == readcorrectedbp[r[3].split('-.-')[1]][r[3].split('-.-')[2].split('/')[0]][0]:
							currdiff = readcorrectedbp[r[3].split('-.-')[1]][r[3].split('-.-')[2].split('/')[0]][1] - int(r[x])
							tempsizes = r[10].split(',')
							tempstarts = r[11].split(',')
							# if c <= 10: print(x, currdiff, tempsizes, tempstarts)
							if x == 1: #left side of gene
								tempsizes[0] = str(int(tempsizes[0]) + currdiff)
								for s in range(1, len(tempstarts)-1):
									tempstarts[s] = str(int(tempstarts[s]) + currdiff)
								r[11] = ','.join(tempstarts)
							elif x == 2: #right side of gene
								tempsizes[-2] = str(int(tempsizes[-2]) + currdiff)
							r[10] = ','.join(tempsizes)
							# if c <= 10: print(tempsizes, tempstarts, r[10:12])
							r[x] = str(readcorrectedbp[r[3].split('-.-')[1]][r[3].split('-.-')[2].split('/')[0]][1])
							r[x+5] = str(readcorrectedbp[r[3].split('-.-')[1]][r[3].split('-.-')[2].split('/')[0]][1])
					reads.write('\t'.join(r) + '\n')
	fusions.close()
	reads.close()
	# import scipy.stats as sps
	# import matplotlib.pyplot as plt
	# # tempa, tempb = max(fastqdistvsbp[0]), max(fastqdistvsbp[1])
	# # fastqdistvsbp[0] = [x/tempa for x in fastqdistvsbp[0]]
	# # fastqdistvsbp[1] = [x / tempb for x in fastqdistvsbp[1]]
	# # print(fastqdistvsbp)
	# print(sps.spearmanr(fastqdistvsbp[0], fastqdistvsbp[1]))
	# fig1 = plt.figure(1)
	# f, ax = plt.subplots()
	# plt.xlabel('fastq distance (norm)')
	# plt.ylabel('dist to breakpoint (norm)')
	# points = ax.scatter([x+1 for x in fastqdistvsbp[0]], [x+1 for x in fastqdistvsbp[1]], alpha=.2)
	# plt.xscale('log')
	# plt.yscale('log')
	# plt.xlim([1,10000])
	# plt.ylim([1,1000000])
	# plt.savefig('fastq-bp-dist-corr.png', dpi=600)




# 'reads': [
# 	['chr12', '51391316', '51489951', 'jointest-SLC4A8-CAPN7-0-0-0--SLC4A8/+/51391317/51515763', '60', '+', '51391316', '51489951', '255,0,0', '19', '172,82,147,136,161,189,92,158,88,147,101,175,134,246,106,162,114,162,252,', '0,49391,59559,60807,62222,66034,67242,68634,69887,70993,72297,78297,79075,79970,83025,83728,94470,97382,98383,'],
# 	['chr12', '51391316', '51489951', 'jointest-SLC4A8-CAPN7-0-0-1--SLC4A8/+/51391317/51515763', '60', '+', '51391316', '51489951', '255,0,0', '19', '172,82,147,136,161,189,92,158,88,147,101,175,134,246,106,162,114,162,252,', '0,49391,59559,60807,62222,66034,67242,68634,69887,70993,72297,78297,79075,79970,83025,83728,94470,97382,98383,'],
# 	['chr12', '51391316', '51489951', 'jointest-SLC4A8-CAPN7-0-0-2--SLC4A8/+/51391317/51515763', '60', '+', '51391316', '51489951', '255,0,0', '19', '172,82,147,136,161,189,92,158,88,147,101,175,134,246,106,162,114,162,252,', '0,49391,59559,60807,62222,66034,67242,68634,69887,70993,72297,78297,79075,79970,83025,83728,94470,97382,98383,'],
# 	['chr12', '51391316', '51489951', 'jointest-SLC4A8-CAPN7-0-0-3--SLC4A8/+/51391317/51515763', '60', '+', '51391316', '51489951', '255,0,0', '19', '172,82,147,136,161,189,92,158,88,147,101,175,134,246,106,162,114,162,252,', '0,49391,59559,60807,62222,66034,67242,68634,69887,70993,72297,78297,79075,79970,83025,83728,94470,97382,98383,'],
# 	['chr12', '51391316', '51489951', 'jointest-SLC4A8-CAPN7-0-0-4--SLC4A8/+/51391317/51515763', '60', '+', '51391316', '51489951', '255,0,0', '19', '172,82,147,136,161,189,92,158,88,147,101,175,134,246,106,162,114,162,252,', '0,49391,59559,60807,62222,66034,67242,68634,69887,70993,72297,78297,79075,79970,83025,83728,94470,97382,98383,']],
# 'CAPN7': {'reads': [
# 	['chr3', '15227836', '15252916', 'jointest-SLC4A8-CAPN7-0-0-0--CAPN7/+/15206152/15252916', '60', '+', '15227836', '15252916', '255,0,0', '15', '129,86,94,147,107,121,145,100,136,76,146,63,131,93,1801,', '0,1137,2605,4682,6030,7188,12636,12917,13616,14341,17689,18895,19490,23094,23279,'],
# 	['chr3', '15227836', '15252916', 'jointest-SLC4A8-CAPN7-0-0-1--CAPN7/+/15206152/15252916', '60', '+', '15227836', '15252916', '255,0,0', '15', '129,86,94,147,107,121,145,100,136,76,146,63,131,93,1801,', '0,1137,2605,4682,6030,7188,12636,12917,13616,14341,17689,18895,19490,23094,23279,'],
# 	['chr3', '15227836', '15252916', 'jointest-SLC4A8-CAPN7-0-0-2--CAPN7/+/15206152/15252916', '60', '+', '15227836', '15252916', '255,0,0', '15', '129,86,94,147,107,121,145,100,136,76,146,63,131,93,1801,', '0,1137,2605,4682,6030,7188,12636,12917,13616,14341,17689,18895,19490,23094,23279,'], ['chr3', '15227836', '15252916', 'jointest-SLC4A8-CAPN7-0-0-3--CAPN7/+/15206152/15252916', '60', '+', '15227836', '15252916', '255,0,0', '15', '129,86,94,147,107,121,145,100,136,76,146,63,131,93,1801,', '0,1137,2605,4682,6030,7188,12636,12917,13616,14341,17689,18895,19490,23094,23279,'], ['chr3', '15227836', '15252916', 'jointest-SLC4A8-CAPN7-0-0-4--CAPN7/+/15206152/15252916', '60', '+', '15227836', '15252916', '255,0,0', '15', '129,86,94,147,107,121,145,100,136,76,146,63,131,93,1801,', '0,1137,2605,4682,6030,7188,12636,12917,13616,14341,17689,18895,19490,23094,23279,']]


# # 	fusionreads = {}
# # 	readToLoc = {}  # We will use this later
# # 	for fusion in new_fusions_found:
# #
# # 	file:
# # 	Line = line.split(‘\t’)
# # 	fusionName = line[3].split(‘-.-’)[0]
# # 	Gene = line[3].split(stuff)
# # 	Readname = line[3].split(stuff)
# # 	fusionreads[fusionName] = {}
# # 	** Option1: fusionreads[fusionname][gene] = [[read1name, read1start, read1end], [read2name, read2start, read2end]]
# # 	** Option2: fusionreads[fusionname][gene] = {leftlocs: [read1start, read2start], rightlocs: [read1end, read2end],
# # 												 readnames: [read1name, read2name]}
# # 	** *Go
# # 	with option2 please
# # 	If
# # 	readname not in readToLoc: readToLoc[readname] = [fusionname]
# # # readToLoc[readname].append([readstart, readend, +/-])
# #
# # For
# # fusion in fusionreads:
# # Closest = []
# # minDistToGeneStart = []
# # For
# # gene in fusionreads[fusion]:
# # Mindist = 10000000000
# # For
# # i in [mean(leftlocs), mean(rightlocs)]:  ##We'll need a much better way to get this than mean but its ok for now
# # For
# # j in geneanno[gene]:  # iterate through transcripts
# # Annostart = j[0] if j[2] == ‘+’ else j[1]
# # If
# # abs(annostart - i) < mindist: mindist = abs(annostart - i)
# # minDistToGeneStart.append([mindist, gene])
# # minDistToGeneStart.sort
# # Transcribed = True if minDistToGeneStart[0][
# # 						  0] < 50 else False  # The 50 is totally arbitrary, we should do tests to figure out a better number
# # print(fusion, transcribed, minDistToGeneStart[0][1])

	# orgFusions = []
	# allMatches = []
	# readNames = {}
	# avgQualScore = 0
	# print('filtering fusions and detecting breakpoints')
	# c = 0
	# #print(new_fusions_found)
	# for i in new_fusions_found:
	# 	supportCount = len(new_fusions_found[i]['readNames'])
	# 	mapScore = round((new_fusions_found[i]['mapScore']/float(supportCount * len(list(i))))/maxMapQ, 3)
	# 	avgBreakpointAgg = 0
	# 	isMito = False
	# 	for loc in new_fusions_found[i]:
	# 		if loc not in ['mapScore', 'readNames']:
	# 			if new_fusions_found[i][loc]['chr'] == "chrM": isMito=True
	# 	if supportCount < int(args.l): metadata.append(['--'.join(list(i)), 'readsup', str(supportCount), new_fusions_found[i]['readNames']])
	# 	# if not ((mapScore >= avgMapQ/float(maxMapQ) or mapScore > 0.8) and mapScore > .5): metadata.append(['--'.join(list(i)), 'mapScore', str(mapScore), new_fusions_found[i]['readNames']])
	# 	if isMito: metadata.append(['--'.join(list(i)), 'isMito', 'M', new_fusions_found[i]['readNames']])
	# 	if (supportCount >= int(args.l) and (args.y or ('chrM' not in i and not isMito))) or i in clinicalF:
	# 		currFusion = [i, str(supportCount), str(mapScore)]
	# 		distTo5 = []
	# 		locInfo = {}
	# 		dupLocusFrac = []
	# 		for loc in new_fusions_found[i]:
	# 			if loc not in ['mapScore', 'readNames', 'repeatScore']:
	# 				leftMed, rightMed = int(median(new_fusions_found[i][loc]['left'])), int(median(new_fusions_found[i][loc]['right']))
	# 				locInfo[loc] = {'chr':new_fusions_found[i][loc]['chr'].strip('chr'), 'loc':leftMed}
	# 				leftFracClose, rightFracClose = 0, 0
	# 				for j in new_fusions_found[i][loc]['left']:
	# 					if abs(j-leftMed) <= 10: leftFracClose += 1
	# 				for j in new_fusions_found[i][loc]['right']:
	# 					if abs(j-rightMed) <= 10: rightFracClose += 1
	# 				leftFracClose = leftFracClose/float(supportCount)
	# 				rightFracClose = rightFracClose/float(supportCount)
	# 				#print(i, loc, leftFracClose, rightFracClose)
	# 				if loc in allGeneLoc:
	# 					if allGeneLoc[loc] != None:
	# 						geneStart = int(allGeneLoc[loc][0][1]) if allGeneLoc[loc][0][3] == '+' else int(allGeneLoc[loc][0][2])
	# 					else: geneStart=None
	# 				else: geneStart=None
	# 				#s = 'm' if Counter(new_fusions_found[i][loc]['strand']).most_common(1)[0][0] == '-' else 'p'
	# 				if loc in geneReads.keys():
	# 					dupLocusFrac.append(round(float(len([a for a in set(geneReads[loc]) if a in new_fusions_found[i]['readNames']])) / len(geneReads[loc]), 3))
	# 				else:
	# 					dupLocusFrac.append('/')
	# 				if leftFracClose > rightFracClose:
	# 					currFusion.append(loc + '-' + new_fusions_found[i][loc]['chr'] + '-' + str(leftMed) + '-' + str(dupLocusFrac[-1]))
	# 					avgBreakpointAgg += leftFracClose
	# 					locInfo[loc]['side'] = 'l'
	# 					if geneStart: distTo5.append(abs(geneStart-rightMed))
	# 					else: distTo5.append(None)
	# 				else:
	# 					currFusion.append(loc + '-' + new_fusions_found[i][loc]['chr'] + '-' + str(rightMed) + '-' + str(dupLocusFrac[-1]))
	# 					avgBreakpointAgg += rightFracClose
	# 					locInfo[loc]['side'] = 'r'
	# 					if geneStart: distTo5.append(abs(geneStart-leftMed))
	# 					else: distTo5.append(None)
	# 		if all(isinstance(x, float) for x in dupLocusFrac):
	# 			currFusion.insert(3, str(round(sum(dupLocusFrac)/len(dupLocusFrac), 3)))
	# 		else: currFusion.insert(3, '/')
	# 		avgBreakpointAgg = avgBreakpointAgg/float(len(new_fusions_found[i].keys()))
	# 		currFusion.append(mapScore*supportCount*avgBreakpointAgg)#*(1-repeatScore))
	# 		avgQualScore += mapScore*supportCount*avgBreakpointAgg#*(1-repeatScore)
	# 		####REMOVING DETERMINATION OF LEFT AND RIGHT GENE
	# 		if None not in distTo5:
	# 			#print(distTo5)
	# 			if distTo5[0] < distTo5[1]:
	# 				temp = currFusion[4]
	# 				currFusion[4] = "3'-" + currFusion[5]
	# 				currFusion[5] = "5'-" + temp
	# 			else:
	# 				currFusion[4] = "3'-" + currFusion[4]
	# 				currFusion[5] = "5'-" + currFusion[5]
	# 				currFusion[-1] *= -1
	# 			if (currFusion[4].split('-')[-3] != currFusion[5].split('-')[-3] or abs(int(currFusion[4].split('-')[-2]) - int(currFusion[5].split('-')[-2])) > args.b or i in clinicalF):
	# 				for j in new_fusions_found[i]['readNames']:
	# 					readNames[j] = {'fusion': i, **copy.deepcopy(locInfo)}
	# 				allMatches += new_fusions_found[i]['readNames']
	# 				orgFusions.append(currFusion)
	# 			else: metadata.append(['--'.join(list(i)), 'tooClose', 'tc', new_fusions_found[i]['readNames']])
	# 		else: metadata.append(['--'.join(list(i)), 'tooClose', 'tc', new_fusions_found[i]['readNames']])
	# correctQ = prefix if os.path.exists(prefix + '.bam') else prefix + '.aligned'
	#
	# if args.s:
	# 	print('making sam file')
	# 	process = subprocess.Popen('samtools view -h -o ' + correctQ + '.sam ' + correctQ + '.bam',stdout=subprocess.PIPE, shell=True)
	# 	print(process.communicate()[0].strip())
	# correctQ = prefix if os.path.exists(prefix + '.sam') else prefix + '.aligned'
	# sam = open(correctQ + '.sam', 'r')
	#
	# allMatches = set(allMatches)
	# readLength = {}
	# sam = open(correctQ + '.sam', 'r')
	# print('checking multi-mapping distance - searching sam file')
	# readSams = {i:{'dist':[],'len':[]} for i in allMatches}
	# c = 0
	# allLenDiff = 0
	# for line in sam:
	# 	line = line.rstrip().split('\t')
	# 	if line[0] in allMatches:
	# 		c += 1
	# 		locs = []
	# 		i = 0
	# 		while line[5][i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']:
	# 			i += 1
	# 		if line[5][i] == 'M': locs.append(0)
	# 		else: locs.append(int(line[5][:i]))
	# 		i = -2
	# 		while line[5][i] not in ['M', 'D', 'I', 'S', 'H', 'X', 'P', 'N']:
	# 			i -= 1
	# 		#TEMPORARY - ASSUMPTIONS
	# 		if line[0] not in readLength:
	# 			readLength[line[0]] = len(line[9])
	# 		#DONE
	# 		if line[5][-1] == 'M': locs.append(readLength[line[0]])
	# 		else: locs.append(readLength[line[0]] - int(line[5][i+1:-1]))
	# 		#New 10/2021
	# 		i = 0
	# 		last = 0
	# 		total = 0
	# 		while i < len(line[5]):
	# 			if not line[5][i].isnumeric():
	# 				if  line[5][i] not in ['S', 'H', 'I', 'P']:
	# 					total += int(line[5][last:i])
	# 				last = i + 1
	# 			i += 1
	# 		fusInfo = readNames[line[0]]
	# 		currDist = 1000000000000000000000
	# 		currGene = None
	# 		for loc in fusInfo.keys():
	# 			if loc != 'fusion':
	# 				if 'chr' + fusInfo[loc]['chr'] == line[2]:
	# 					temp = abs(fusInfo[loc]['loc'] - int(line[3]))
	# 					if temp < currDist:
	# 						currGene = loc
	# 						currDist = temp
	# 		lenDiff=1
	# 		if currGene != None:
	# 			if currGene[:3] != 'chr' and currGene in allGeneLoc:
	# 				lenDiff = total/(allGeneLoc[currGene][0][2]-allGeneLoc[currGene][0][1])
	# 		allLenDiff += lenDiff
	# 		#DONE
	# 		if line[1] == '16' or line[1] == '2064':
	# 			locs[0] = readLength[line[0]]-locs[0]
	# 			locs[1] = readLength[line[0]] - locs[1]
	# 			locs = locs[::-1]
	# 		locs.append(line[2] + '-' + line[3])
	# 		readSams[line[0]]['dist'].append(locs)
	# 		readSams[line[0]]['len'].append(lenDiff)
	# sam.close()
	# c = 0
	# fastqFusionLocs = {}
	# keysToRemove = []
	#
	# for i in readSams:
	# 	c += 1
	# 	if len(readSams[i]['dist']) == 2:
	# 		readSams[i]['dist'].sort()
	# 		readSams[i]['dist'].append(readSams[i]['dist'][1][0]-readSams[i]['dist'][0][1])
	# 		if readNames[i]['fusion'] not in fastqFusionLocs:
	# 			fastqFusionLocs[readNames[i]['fusion']] = {'loc':[], 'len':[], 'fastqCov':[]}
	# 		fastqFusionLocs[readNames[i]['fusion']]['loc'].append(readSams[i]['dist'])
	# 		tempavg = sum(readSams[i]['len'])/float(len(readSams[i]['len']))
	# 		fastqFusionLocs[readNames[i]['fusion']]['len'].append(tempavg + ((max(readSams[i]['len'])-tempavg)/2))#sum(readSams[i]['len'])/float(len(readSams[i]['len'])))
	# 		fastqFusionLocs[readNames[i]['fusion']]['fastqCov'].append(((readSams[i]['dist'][1][1]-readSams[i]['dist'][0][0])-(readSams[i]['dist'][1][0]-readSams[i]['dist'][0][1]))/readLength[i])
	# 	#THIS IS TEMPORARY< REMOVE FOR MORE FLEXIBILITY/3-GENE
	# 	else: keysToRemove.append(i)
	# 	#DONE
	#
	# for i in fastqFusionLocs:
	# 	temp = [a[-1] for a in fastqFusionLocs[i]['loc']]
	# 	temp.sort()
	# 	fastqFusionLocs[i]['loc'].append(temp[int(len(temp) / 2)])
	#
	#
	# # print('fusions filtered')
	# reads = open(outfilename + "Reads.bed", "w")
	# fusions = open(outfilename + "Fusions.tsv", "w")
	# if len(orgFusions) > 0:
	# 	avgQualScore = avgQualScore/len(orgFusions)
	# orgFusions.sort(key=lambda x:abs(x[-1]), reverse=True)
	# fusions.write("#name\tspanning reads\tmapping score(1 is good)\tavg frac of reads at loci in fusion\t3' breakpoint\t5' breakpoint\n")
	# fusionReadLocs = {}
	# orgFusionsDict = {}
	#
	# for i in orgFusions:
	# 	if abs(i[-1]) > avgQualScore*.01:
	# 		orgFusionsDict[i[0]] = i[:-1]
	# 		fusionReadLocs[i[0]] = {i[5].split('-')[1] + '-' + i[5].split('-')[-3]:[], i[4].split('-')[1] + '-' + i[4].split('-')[-3]:[]}
	#
	# count = 0
	# last = ""
	# print("num pre final filter", len(fusionReadLocs.keys()))
	# bedLinesFiltered = []
	# fusionDist = {a:[] for a in list(fusionReadLocs.keys())}
	# for line in bedLines:
	# 	thisName = line[3].split('--')[0]
	# 	if thisName in allMatches:
	# 		fName = readNames[thisName]['fusion'] ##Maybe this is where we figure out the orientation of the fusion name
	# 		currName = line[3].split('--')[1].split('/')[0]
	# 		if currName[:3] == 'chr':
	# 			#currName = currName.split('.')[0]
	# 			for loc in readNames[thisName]:
	# 				if loc[:3] == 'chr':
	# 					if abs(float(currName.split('-')[1]) - float(loc.split('-')[1])) <= args.b * 3:
	# 						currName = loc
	# 		locName = line[3].split('--')[1].split('/')[0] if '-' not in line[3].split('--')[1].split('/')[0] else line[3].split('--')[1].split('/')[0].split('.')[0]
	# 		if fName in fusionReadLocs:
	# 			if locName + '-' + line[0] in fusionReadLocs[fName]:
	# 				fusionReadLocs[fName][locName + '-' + line[0]].append([int(line[1]), int(line[2])])
	# 		if currName in readNames[thisName].keys():
	# 			###DEFINITELY NEED ORIENTATION OF GENES IN FUSION BY HERE
	# 			line[3] = '-.-'.join(['--'.join(list(fName))] + line[3].split('--'))
	# 			if line[-2][-1] != ',':
	# 				line[-2] += ','
	# 				line[-1] += ','
	# 			bedLinesFiltered.append(line)
	# c = 0
	# x = 0
	# finalFusions = []
	# almostDone = []
	# adGenes = []
	# for fusion in fusionReadLocs:
	# 	fusionlist = list(fusion)
	# 	wasWritten = False
	# 	theseLocs = []
	# 	SSdist = []
	# 	readEnds = []
	# 	counts = []
	# 	chrs = []
	# 	for loc in fusionReadLocs[fusion]:
	# 		chrs.append(loc.split('-')[-1])
	# 		temp = [list(i) for i in zip(*fusionReadLocs[fusion][loc])]
	# 		far = True
	# 		if len(temp) > 0:
	# 			temp[0].sort()
	# 			temp[1].sort()
	# 			counts.append(min(temp[0]))
	# 			counts.append(max(temp[1]))
	# 			leftSS1 = binarySearch(junctions[loc.split('-')[-1]], temp[0][int(len(temp[0])/2)])
	# 			rightSS1 = binarySearch(junctions[loc.split('-')[-1]], temp[1][int(len(temp[1]) / 2)])
	# 			if len(temp[0]) % 2 == 0 and len(temp) > 2:
	# 				leftSS1 = binarySearch(junctions[loc.split('-')[-1]], temp[0][int(len(temp[0]) / 2)]+1)
	# 				rightSS1 = binarySearch(junctions[loc.split('-')[-1]], temp[1][int(len(temp[1]) / 2)]+1)
	# 				leftDistToSS = min([abs(leftSS1-temp[0][int(len(temp[0])/2)]), abs(leftSS2-temp[0][int(len(temp[0])/2)+1])])
	# 				rightDistToSS = min([abs(rightSS1 - temp[1][int(len(temp[1]) / 2)]), abs(rightSS2 - temp[1][int(len(temp[1]) / 2) + 1])])
	# 			else:
	# 				leftDistToSS = abs(leftSS1-temp[0][int(len(temp[0])/2)])
	# 				rightDistToSS = abs(rightSS1-temp[1][int(len(temp[1])/2)])
	# 			SSdist.append(sorted([[leftDistToSS, temp[0][int(len(temp[0])/2)]], [rightDistToSS,temp[1][int(len(temp[1])/2)]]]))
	# 			theseLocs.append(loc.split('-')[0])
	# 		else: far = False
	# 	if len(counts) > 2:
	# 		if abs(counts[0] - counts[2]) < args.b or abs(counts[0] - counts[3]) < args.b or \
	# 				abs(counts[1] - counts[2]) < args.b or abs(counts[1] - counts[3]) < args.b:
	# 			far = False
	# 	fastqDist = 500
	#
	# 	if fusion in fastqFusionLocs:
	# 		fastqDist = abs(fastqFusionLocs[fusion]['loc'][-1])
	# 		lenDiff = str(sum(fastqFusionLocs[fusion]['len']) / len(fastqFusionLocs[fusion]['len']))
	# 		fastqCov = str(sum(fastqFusionLocs[fusion]['fastqCov']) / len(fastqFusionLocs[fusion]['fastqCov']))
	# 	if fastqDist >= 15: metadata.append(['--'.join(fusionlist), 'fastqDist', str(fastqDist), new_fusions_found[fusion]['readNames']])#metadata["fastqDist"] += 1
	# 	# print(list(set(chrs)))
	# 	if not (len(list(set(chrs))) > 1 or far): metadata.append(['--'.join(fusionlist), 'tooClose2', 'tc2', new_fusions_found[fusion]['readNames']])#metadata["tooClose2"] += 1
	# 	if (fastqDist < 15 and (len(list(set(chrs))) > 1 or far)) or fusion in clinicalF:
	# 		if len(SSdist) > 1:
	# 			if SSdist[0][0][0] <= 10 and SSdist[1][0][0] <= 10:
	# 				if theseLocs[0] in allGeneLoc:
	# 					if allGeneLoc[theseLocs[0]] != None:
	# 						aDistToProm = abs(SSdist[0][1][1] - int(allGeneLoc[theseLocs[0]][0][1])) if allGeneLoc[theseLocs[0]][0][3] == '+' else abs(SSdist[0][1][1] - int(allGeneLoc[theseLocs[0]][0][2]))
	# 					else:
	# 						aDistToProm = 1000000000000000
	# 				else: aDistToProm = 1000000000000000
	# 				if theseLocs[1] in allGeneLoc:
	# 					if allGeneLoc[theseLocs[1]] != None:
	# 						bDistToProm = abs(SSdist[1][1][1] - int(allGeneLoc[theseLocs[1]][0][1])) if allGeneLoc[theseLocs[1]][0][3] == '+' else abs(SSdist[0][1][1] - int(allGeneLoc[theseLocs[1]][0][2]))
	# 					else:
	# 						bDistToProm = 1000000000000000
	# 				else: bDistToProm = 1000000000000000
	# 				temp = [a.split('-') for a in list(orgFusionsDict[fusion][-2:])]
	# 				if temp[0][-4] == theseLocs[0]:
	# 					temp[0][-2] = str(SSdist[0][0][1])
	# 					temp[1][-2] = str(SSdist[1][0][1])
	# 				elif temp[1][-4] == theseLocs[0]:
	# 					temp[0][-2] = str(SSdist[1][0][1])
	# 					temp[1][-2] = str(SSdist[0][0][1])
	# 				orgFusionsDict[-2] = '-'.join(temp[0])
	# 				orgFusionsDict[-1] = '-'.join(temp[1])
	# 				orgFusionsDict[fusion].append(lenDiff)
	# 				orgFusionsDict[fusion].append(fastqCov)
	# 				almostDone.append(orgFusionsDict[fusion])
	# 				adGenes.append(fusionlist[0])
	# 				adGenes.append(fusionlist[1])
	# 				c += 1
	# 			else: metadata.append(['--'.join(fusionlist), 'ssDist', '-'.join([str(x) for x in SSdist]), new_fusions_found[fusion]['readNames']])#metadata["ssdist"] += 1
	# 		elif len(SSdist) == 1:
	# 			if SSdist[0][0][0] < 10:
	# 				c += 1
	# 				almostDone.append(orgFusionsDict[fusion])
	# 				adGenes.append(fusionlist[0])
	# 				adGenes.append(fusionlist[1])
	# 			else: metadata.append(['--'.join(fusionlist), 'ssDist', '-'.join([str(x) for x in SSdist]), new_fusions_found[fusion]['readNames']])#metadata["ssdist"] += 1
	# 		else: metadata.append(['--'.join(fusionlist), 'ssDist', '-'.join([str(x) for x in SSdist]), new_fusions_found[fusion]['readNames']])#metadata["ssdist"] += 1
	# c = 0
	# freq = Counter(adGenes)
	# for temp in almostDone:
	# 	good = True
	# 	fusion, fusionlist = temp[0], list(temp[0])
	# 	if len(temp) > 8 and temp[0] not in clinicalF:
	# 		if float(temp[-1]) < 0.7:
	# 			good = False
	# 			metadata.append(['--'.join(fusionlist), 'fastqCov', str(temp[-1]), new_fusions_found[fusion]['readNames']])
	# 		if float(temp[-2]) > 0.9:
	# 			good = False
	# 			metadata.append(['--'.join(fusionlist), 'geneCov', str(temp[-2]), new_fusions_found[fusion]['readNames']])
	# 		if args.w:
	# 			if freq[fusionlist[0]] > 2 or freq[fusionlist[1]] > 2:
	# 				good = False
	# 				metadata.append(['--'.join(fusionlist), 'repeat', 'r', new_fusions_found[fusion]['readNames']])
	# 	if good:
	# 		finalFusions.append(temp[0])
	# 		##NEED TO ADD FUSION GENE ORDER HERE
	# 		fusions.write('\t'.join(['--'.join(list(temp[0]))] + temp[1:7]) + '\n')
	# 		c += 1
	# print("final output", c)
	# fusions.close()
	# print('fusions written')
	# for i in metadata:
	# 	i[3] = ','.join(i[3])
	# 	meta.write('\t'.join(i) + '\n')#i + "\t" + str(metadata[i]) + "\n")
	# meta.close()
	# finalFusions = set(finalFusions)
	# printNames = open(outfilename + "readNames.txt", "w")
	# for line in bedLinesFiltered:
	# 	if line[3].split('-.-')[0] in finalFusions:
	# 		reads.write("\t".join(line) + "\n")
	# 		printNames.write(line[3].split('-.-')[1] + '\n')
	# reads.close()
	# printNames.close()
	# print('reads written')
	#


	# subprocess.run(['perl', os.path.dirname(os.path.realpath(__file__)) + '/bed12ToGTF.pl'],
	# 			   stdin=open(args.o + prefix + 'Reads.bed', 'r'), stdout=open(args.o + prefix + 'Reads.gtf', 'w'))
# if int(args.k) > 0:
# 	readNames = []
# 	print('remapping reads to fusions')
# 	with open(outfilename + "readNames.txt", 'r') as names:
# 		for line in names:
# 			readNames.append(line.rstrip())
# 	readNames = set(readNames)
# 	#Filter reads to only the double mapped reads
# 	leadingChar = '@' if args.r.split('.')[-1] == 'fastq' or args.r.split('.')[-1] == 'fq' else '>'
# 	with open(args.r, 'r') as reads, open(outfilename + "Filtered.fa", "w") as faOut:
# 		writeRead = False
# 		c, d = 0, 0
# 		for line in reads:
# 			c += 1
# 			if c % 20000000 == 0: print(c, d)
# 			if line[0] == leadingChar:
# 				if line.rstrip('\n').lstrip(leadingChar) in readNames:
# 					writeRead = True
# 					d += 1
# 					if leadingChar == '@': faOut.write(">" + line.lstrip('@'))
# 					else: faOut.write(line)
# 				else:
# 					writeRead = False
# 			elif writeRead:
# 				faOut.write(line)
# 				writeRead = False
# 	print('reads filtered')
# 	fusions = {}
# 	firstLine = []
# 	with open(outfilename + "Fusions.tsv", 'r') as thesefusions, open(outfilename + 'Locs.bed', 'w') as bedFile:
# 		for line in thesefusions:
# 			line = line.rstrip().split('\t')
# 			if line[0][0] != '#' and len(line) > 6:
# 				chr1, center1 = line[5].split('-')[-3:-1]
# 				chr2, center2 = line[6].split('-')[-3:-1]
# 				if int(center1) > args.k and int(center2) > args.k:
# 					name1 = line[0] + '->' + '-'.join(line[5].lstrip("3'-").split('-')[:-2])
# 					name2 = line[0] + '->' + '-'.join(line[6].lstrip("5'-").split('-')[:-2])
# 					fusions[line[0]] = {'line':line}
# 					fusions[line[0]]['-'.join(line[5].lstrip("3'-").split('-')[:-2])] = {'side':"3'", 'chr':chr1, 'bp':int(center1), 'left':[], 'right':[], 'reads':[], 'mapQ':[]}
# 					fusions[line[0]]['-'.join(line[6].lstrip("5'-").split('-')[:-2])] = {'side':"5'", 'chr':chr2, 'bp':int(center2), 'left':[], 'right':[], 'reads':[], 'mapQ':[]}
# 					bedFile.write('\t'.join([chr1, str(int(center1)-args.k), str(int(center1) + args.k), name1]) + '\n')
# 					bedFile.write('\t'.join([chr2, str(int(center2)-args.k), str(int(center2) + args.k), name2]) + '\n')
# 			elif line[0][0] == '#':
# 				firstLine = line
# 	process = subprocess.Popen('bedtools getfasta -fi ' + args.g + ' -bed ' + outfilename + 'Locs.bed' + ' -fo ' + outfilename + 'Genome.fa' + ' -name; ' +
# 							   'minimap2 -a ' + outfilename + 'Genome.fa ' + outfilename + "Filtered.fa" + ' > ' + outfilename + 'Remapped.sam; ' +
# 							   "sam2bed < " + outfilename + 'Remapped.sam' + ' > ' + outfilename + 'Remapped-unfilt.bed',stdout=subprocess.PIPE, shell=True)
# 	print(process.communicate()[0].strip())
# 	maxMapQ = 0
# 	with open(outfilename + 'Remapped-unfilt.bed', 'r') as remapped:
# 		for line in remapped:
# 			line = line.rstrip().split('\t')
# 			if int(line[4]) > maxMapQ: maxMapQ = int(line[4])
# 			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['mapQ'].append(int(line[4]))
# 			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['reads'].append(line[3])
# 			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['left'].append(int(line[1])-args.s)
# 			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['right'].append(int(line[2])-args.s)
# 	newFusions = open(outfilename + 'FusionsRemapped.tsv', 'w')
# 	firstLine.insert(2, 'confirmed reads')
# 	newFusions.write('\t'.join(firstLine) + '\n')
# 	doubleMappedReads = []
# 	for fusion in fusions:
# 		good = 0
# 		fusionReads = []
# 		mapScores = []
# 		for loc in fusions[fusion]:
# 			if loc not in ['line', 'conf reads']:
# 				if len(fusions[fusion][loc]['left']) > 0:
# 					leftAvg, rightAvg = sum(fusions[fusion][loc]['left'])/len(fusions[fusion][loc]['left']), \
# 										sum(fusions[fusion][loc]['right'])/len(fusions[fusion][loc]['right'])
# 					if abs(leftAvg + rightAvg) > 50:
# 						good += 1
# 						# if abs(leftAvg) > abs(rightAvg):
# 						# 	fusions[fusion][loc]['bp'] += max(fusions[fusion][loc]['right'])
# 						# else:
# 						# 	fusions[fusion][loc]['bp'] += min(fusions[fusion][loc]['left'])
# 						fusionReads.append(fusions[fusion][loc]['reads'])
# 						mapScores += fusions[fusion][loc]['mapQ']
# 						# if fusions[fusion][loc]['side'] == "3'":
# 						# 	fusions[fusion]['line'][5] = '-'.join(fusions[fusion]['line'][5].split('-')[:-1] + [str(fusions[fusion][loc]['bp'])])
# 						# else:
# 						# 	fusions[fusion]['line'][6] = '-'.join(fusions[fusion]['line'][6].split('-')[:-1] + [str(fusions[fusion][loc]['bp'])])
# 		if good >= 2:
# 			shared = list(set(fusionReads[0]) & set(fusionReads[1]))
# 			if len(shared) > 0:#good >= 2:
# 				for i in shared: doubleMappedReads.append(fusion + i)
# 				fusions[fusion]['line'].insert(2, str(len(shared)))
# 				fusions[fusion]['line'][3] = str(round((sum(mapScores)/len(mapScores))/float(maxMapQ), 3))
# 				newFusions.write('\t'.join(fusions[fusion]['line']) + '\n')
# 	newFusions.close()
# 	doubleMappedReads = set(doubleMappedReads)
# 	goodReadNames = []
# 	with open(outfilename + 'Remapped-unfilt.bed', 'r') as remapped, open(outfilename + 'Remapped.bed', 'w') as remapFilt:
# 		for line in remapped:
# 			line = line.split('\t')
# 			if line[0].split('->')[0] + line[3] in doubleMappedReads:
# 				remapFilt.write('\t'.join(line))
# 				goodReadNames.append(line[3])
# 	readsOut = open(outfilename + 'RemappedReads.bed', 'w')
# 	goodReadNames = set(goodReadNames)
# 	for line in open(outfilename + "Reads.bed", 'r'):
# 		temp = line.split('\t')
# 		if temp[3].split('-.-')[1] in goodReadNames:
# 			readsOut.write(line)
#
#
# 	process = subprocess.Popen('python ' + os.path.dirname(os.path.realpath(__file__)) + '/makeAlnSeq.py -f ' + outfilename +
# 							   'Genome.fa -r ' + outfilename + 'Remapped-unfilt.bed; rm ' + outfilename + 'Remapped-unfilt.bed',stdout=subprocess.PIPE, shell=True)
# 	print(process.communicate()[0].strip())


if not args.i:
	# for line in open(outfilename + 'Fusions.tsv', 'r'):
	# 	if line[0] != '#':
	# 		line = line.rstrip().split('\t')
	# 		fusion = line[0].split('--')#sorted(line[0].split('--'))
	# 		# print(fusion)
	# 		left.append(fusion[0])
	# 		right.append(fusion[1])
	# left = set(left)
	# right = set(right)
	# print("left")
	# print(left)
	leftOut = open(outfilename + 'Reads-l.bed', 'w')
	rightOut = open(outfilename + 'Reads-r.bed', 'w')
	leftReads, rightReads = [], []
	c = 0
	freads = []
	for line in open(outfilename + 'Reads.bed', 'r'):
		if line[:5] != 'track':
			line2 = line.strip().split('\t')
			fusionname = line2[3].split('-.-')[0].split('--')
			geneName = line2[3].split('-.-')[2].split('/')[0]
			if geneName[:3] == 'chr': geneName = geneName.split('.')[0]
			if fusionname.index(geneName) == 0: leftOut.write(line)
			else: rightOut.write(line)
			freads.append(line2[3].split('-.-')[1])
		# if geneName in left:
		# 	c += 1
		# 	leftOut.write(line)
		# elif geneName in right:
		# 	rightOut.write(line)
	freads = set(freads)
	leftOut.close()
	rightOut.close()
	freadsout = open(outfilename + '-fusionreads.fastq', 'w')
	last = None
	for line in open(args.r):
		if line[0] == '@' or line[0] == '>': last = line.split(' ')[0].lstrip('@>')
		if last in freads: freadsout.write(line)
	freadsout.close()
	# print(c)
	#print('python3 ' + args.f + ' collapse -g ' + args.g + ' -r ' + args.r + ' --generate_map -q ' + outfilename + 'Reads.bed -o ' + prefix + '.fusions',)
	#print('python3 ' + args.f + ' collapse -g ' + args.g + ' -r ' + args.r + ' --generate_map -q ' + outfilename + 'Reads-l.bed -o ' + prefix + '.fusions.l')
	process = subprocess.Popen(
		#collapse breaks.simplen/he_v2.3.5_pass.fastq --generate_map -q 17-03-2021he_v2.3.5_passReads-1.bed -o he.fusions.collapse
		'python3 ' + args.f + ' collapse --stringent --temp_dir /scratch/cafelton/ -g ' + args.g + ' -r ' + outfilename + '-fusionreads.fastq' + ' --generate_map -q ' + outfilename + 'Reads-l.bed -o ' + prefix + '.fusions.l' +
		'; python3 ' + args.f + ' collapse --stringent --temp_dir /scratch/cafelton/ -g ' + args.g + ' -r ' + outfilename + '-fusionreads.fastq' + ' --generate_map -q ' + outfilename + 'Reads-r.bed -o ' + prefix + '.fusions.r',
		stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())

if not args.j and os.path.exists(prefix + '.fusions.l.isoform.read.map.txt'):
	from fusionalignmentvisualization import *
	readisos = {}
	for line in open(outfilename + 'Reads.bed', 'r'):
		if line[:5] != 'track':
			line = line.split('\t')
			ids = line[3].split('-.-')
			readisos[ids[1]] = [ids[0]]
	for line in open(prefix + '.fusions.l.isoform.read.map.txt', 'r'):
		line = line.rstrip().split('\t')
		line[1] += ',' + line[0].split('-.-')[1]
		for read in line[1].split(','):
			readisos[read].append(line[0])
	for line in open(prefix + '.fusions.r.isoform.read.map.txt', 'r'):
		line = line.rstrip().split('\t')
		line[1] += ',' + line[0].split('-.-')[1]
		for read in line[1].split(','):
			readisos[read].append(line[0])
	fusionisos = {}
	for r in readisos:
		if len(readisos[r]) > 1:
			if readisos[r][0] not in fusionisos: fusionisos[readisos[r][0]] = {}
			if frozenset(readisos[r][1:]) not in fusionisos[readisos[r][0]]: fusionisos[readisos[r][0]][frozenset(readisos[r][1:])] = 0
			fusionisos[readisos[r][0]][frozenset(readisos[r][1:])] += 1
	goodisos, isocount = {}, 0
	fusionout = open(outfilename + 'IsoformFusions.tsv', 'w')
	for line in open(outfilename + 'Fusions.tsv', 'r'):
		line = line.split('\t')
		if line[0] in fusionisos:
			d,s,t = 0,0,0
			for i in fusionisos[line[0]]:
				if len(i) < 2:
					c = 0
					for j in fusionisos[line[0]]:
						if list(i)[0] in j: c += 1
					if fusionisos[line[0]][i] < 3: fusionisos[line[0]][i] = 'lowSupportSingle' + str(fusionisos[line[0]][i])
					elif c > 1: fusionisos[line[0]][i] = 'singleIsoSubsumedInDouble' + str(fusionisos[line[0]][i])
					elif c <= 1 and fusionisos[line[0]][i] >= 3:
						s += 1
						t += fusionisos[line[0]][i]
						if list(i)[0] not in goodisos: goodisos[list(i)[0]] = []
						goodisos[list(i)[0]].append([str(fusionisos[line[0]][i]), 'combIso' + str(isocount), 'single'])
						fusionisos[line[0]][i] = 'combIso' + str(isocount)
						isocount += 1
				elif fusionisos[line[0]][i] >= 3:
					d += 1
					t += fusionisos[line[0]][i]
					for l in list(i):
						if l not in goodisos: goodisos[l] = []
						goodisos[l].append([str(fusionisos[line[0]][i]), 'combIso' + str(isocount), 'double'])
					fusionisos[line[0]][i] = 'combIso' + str(isocount)
					isocount += 1
				else: fusionisos[line[0]][i] = 'lowSupportDouble' + str(fusionisos[line[0]][i])
			line[1] = '/'.join([str(d), str(s), str(t), line[1]])
			fusionout.write('\t'.join(line))
	out = open(outfilename + "ReadsIsoSupport.bed", 'w')
	toColor = {'lowSupportDouble':'109,17,189', 'combIso':'26,189,17', 'singleIsoSubsumedInDouble':'17,166,189', 'lowSupportSingle':'169,17,189', 'None':'189,17,49'}
	readCounts = {'lowSupportDouble':0, 'combIso':0, 'singleIsoSubsumedInDouble':0,'lowSupportSingle':0, 'None':0}
	totreads = 0
	for line in open(outfilename + 'Reads.bed', 'r'):
		if line[:5] != 'track':
			line = line.split('\t')
			ids = line[3].split('-.-')
			fusionToSupport = 'None' if len(readisos[ids[1]]) == 1 else str(fusionisos[ids[0]][frozenset(readisos[ids[1]][1:])])
			line[3] = '-.-'.join([fusionToSupport] + ids)
			line[8] = toColor[fusionToSupport.rstrip('1234567890')]
			readCounts[fusionToSupport.rstrip('1234567890')] += 1
			totreads += 1
			out.write('\t'.join(line))
	for i in readCounts: readCounts[i] = readCounts[i]/totreads
	print(readCounts)
	out.close()
	# # for l in goodisos:
	# # 	print(l)
	# # 	for m in goodisos[l]: print(m)
	fusionout.close()
	readsout = open(outfilename + 'IsoformReads.bed', 'w')
	for file in [prefix + '.fusions.l.isoforms.bed', prefix + '.fusions.r.isoforms.bed']:
		for line in open(file):
			line = line.split('\t')
			if line[3] in goodisos:
				if len(goodisos[line[3]]) > 1:
					for x in range(len(goodisos[line[3]])):
						temp = list(set([y[1] for y in goodisos[line[3]]])-set([goodisos[line[3]][x][1]]))
						temp.sort()
						goodisos[line[3]][x].append("*this-loc-shared-with-" + ','.join(temp))
				for m in goodisos[line[3]]:
					line[3] = '-.-'.join([line[3].split('-.-')[0], '|'.join(m), line[3].split('-.-')[2]])
					readsout.write('\t'.join(line))
	readsout.close()
	# print('rm -f ' + prefix + '.fusions.l.isoforms.fa ' + prefix + '.fusions.l.isoforms.bed ' + prefix + '.fusions.l.isoform.read.map.txt ' + prefix + '.fusions.r.isoforms.fa ' + prefix + '.fusions.r.isoforms.bed ' + prefix + '.fusions.r.isoform.read.map.txt ' + outfilename + 'Reads-l.bed ' + outfilename + 'Reads-r.bed ' + outfilename + '-fusionreads.fastq')
	process = subprocess.Popen('rm -f ' + prefix + '.fusions.l.isoforms.fa ' + prefix + '.fusions.l.isoforms.bed ' + prefix + '.fusions.l.isoform.read.map.txt ' + prefix + '.fusions.r.isoforms.fa ' + prefix + '.fusions.r.isoforms.bed ' + prefix + '.fusions.r.isoform.read.map.txt ' + outfilename + 'Reads-l.bed ' + outfilename + 'Reads-r.bed ' + outfilename + '-fusionreads.fastq',
		stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())
	make_fusion_plots(outfilename + 'IsoformReads.bed', args.t)




	# fusionisos = {}
	# fusionHeader = ''
	# for line in open(outfilename + 'Fusions.tsv', 'r'):
	# 	if line[0] == '#': fusionHeader = line
	# 	else:
	# 		line = line.strip().split('\t')
	# 		# temp = line[0].split('--')
	# 		fusionisos[line[0]] = {'reads':[], 'isos':[], 'line':line}
	# for line in open(outfilename + 'Reads.bed', 'r'):
	# 	if line[:5] != 'track':
	# 		line = line.split('\t')
	# 		info = line[3].split('-.-')
	# 		if info[0] in fusionisos.keys():
	# 			fusionisos[info[0]]['reads'].append(info[1])
	# lIsoSup, rIsoSup = {}, {}
	# lIsoKey, rIsoKey = {}, {}
	# c = 0
	# for line in open(prefix + '.fusions.l.isoform.read.map.txt', 'r'):
	# 	temp = line.strip().split('\t')
	# 	lIsoKey['lIso' + str(c)] = temp[0]
	# 	c += 1
	# 	for i in temp[1].split(','):
	# 		lIsoSup[i] = temp[0]#'lIso' + str(c) #temp[0]
	# c = 0
	# for line in open(prefix + '.fusions.r.isoform.read.map.txt', 'r'):
	# 	temp = line.strip().split('\t')
	# 	rIsoKey['rIso' + str(c)] = temp[0]
	# 	c += 1
	# 	for i in temp[1].split(','):
	# 		rIsoSup[i] = temp[0]#'rIso' + str(c) #temp[0]
	# lSet, rSet = set(lIsoSup.keys()), set(rIsoSup.keys())
	# for i in fusionisos.keys():
	# 	fusionisos[i]['reads'] = set(fusionisos[i]['reads'])
	# 	for j in fusionisos[i]['reads']:
	# 		inl, inr = '.', '.'
	# 		# TESTING
	# 		if j in lSet:
	# 			if lIsoSup[j].split('-.-')[0]==i:
	# 				inl = lIsoSup[j]
	# 		if j in rSet:
	# 			if rIsoSup[j].split('-.-')[0]==i:
	# 				inr = rIsoSup[j]
	# 		if inl != '.' or inr != '.':
	# 			fusionisos[i]['isos'].append('-=-'.join([inl, inr]))
	# 		# END TESTING
	# c = 0
	# #print(lSet)
	# isoToCombName = {}
	# fusionsOut = open(outfilename + 'IsoformFusions.tsv', 'w')
	# fusionsOut.write(fusionHeader)
	# #print(fusionisos)
	# for i in fusionisos.keys():
	# 	combIsoCounts = Counter(fusionisos[i]['isos'])
	# 	#if c < 30: print(i, combIsoCounts)
	# 	totIsos, doubleIsos, readsSupIsos = 0, 0, 0
	# 	for j in combIsoCounts.keys():
	# 		temp = j.split('-=-')
	# 		numsides = 'single'
	# 		if temp[0] != '.' and temp[1] != '.':
	# 			doubleIsos += 1
	# 			numsides = 'double'
	# 		totIsos += 1
	# 		#print(tempr)
	# 		readsSupIsos += combIsoCounts[j]
	# 		if temp[0] != '.':
	# 			if temp[0].split('-.-')[1] not in isoToCombName: isoToCombName[temp[0].split('-.-')[1]] = []
	# 			isoToCombName[temp[0].split('-.-')[1]].append("#" + i)
	# 			isoToCombName[temp[0].split('-.-')[1]].append(str(combIsoCounts[j]) + '|combIso' + str(c) + '|' + numsides)
	# 		if temp[1] != '.':
	# 			if temp[1].split('-.-')[1] not in isoToCombName: isoToCombName[temp[1].split('-.-')[1]] = []
	# 			isoToCombName[temp[1].split('-.-')[1]].append("#" + i)
	# 			isoToCombName[temp[1].split('-.-')[1]].append(str(combIsoCounts[j]) + '|combIso' + str(c) + '|' + numsides)
	# 		c += 1
	# 	fusionisos[i]['line'][1] = str(doubleIsos) + '/' + str(totIsos) + '/' + str(readsSupIsos) + '/' + fusionisos[i]['line'][1]
	# 	if totIsos > 0:
	# 		fusionsOut.write('\t'.join(fusionisos[i]['line']) + '\n')
	# readsOut = open(outfilename + 'IsoformReads.bed', 'w')
	# #TESTING
	# # multiGeneIsos = []
	# # for i in isoToCombName.keys():
	# # 	fusionList = []
	# # 	for j in isoToCombName[i]:
	# # 		if j[0] == '#': fusionList.append(j)
	# # 	if len(Counter(fusionList).keys()) > 1:
	# # 		multiGeneIsos.append([i, fusionList])
	# # print(multiGeneIsos)
	# # print(len(multiGeneIsos))
	# c = 0
	# written = []
	# for i in ['l', 'r']:
	# 	for line in open(prefix + '.fusions.' + i + '.isoforms.bed', 'r'):
	# 		line = line.strip().split('\t')
	# 		info = line[3].split('-.-')
	# 		if info[1] in isoToCombName.keys():
	# 			# c += 1
	# 			# if c < 10:
	# 			# 	print(isoToCombName[info[1]])
	# 			for j in isoToCombName[info[1]]:
	# 				if j[0] != '#':
	# 					temp = info
	# 					info[1] = j
	# 					line[3] = '-.-'.join(info)
	# 					if line[3] not in written:
	# 						written.append(line[3])
	# 						readsOut.write('\t'.join(line) + '\n')
	# fusionsOut.close()
	# readsOut.close()
	# process = subprocess.Popen('rm -f ' + prefix + '.fusions.l.isoforms.fa ' + prefix + '.fusions.l.isoforms.bed ' + prefix + '.fusions.l.isoform.read.map.txt' + prefix + '.fusions.r.isoforms.fa ' + prefix + '.fusions.r.isoforms.bed ' + prefix + '.fusions.r.isoform.read.map.txt ' + outfilename + 'Reads-l.bed ' + outfilename + 'Reads-r.bed',stdout=subprocess.PIPE, shell=True)
	# print(process.communicate()[0].strip())


# read on one side can be matched to multiple reads on the other (in combIsoCounts)
# must pick best
# why are reads from diff. fusions getting intp fusion?
# make isoform support stats lign up
# when multiple fusions in the same gene, it messes up
# check for any fusions with 2 unique loci - see if they look the same on the og version of testing section and new version
