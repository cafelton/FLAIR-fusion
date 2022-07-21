import sys, csv, subprocess, os, argparse
from datetime import date
from statistics import median
from collections import Counter
import copy
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

#Removing all repetitive region references

parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python3 19-03-2021-flair-to-fusions-pipe.py -g genome.fa -t anno.gtf -a anno-short.gtf -f path-to-flair -r reads.fastq')
parser.add_argument('-o', '--output', action='store', dest='o', default=date.today().strftime("%d-%m-%Y"), help='output file name base (default: date)')
parser.add_argument('-r', '--reads', action='store', dest='r', default="", help='.fa or fq file')
parser.add_argument('-m', '--bedFile', action='store', dest='m', default="", help='.bed file')
parser.add_argument('-f', '--flair', action='store', dest='f', default=os.path.dirname(os.path.realpath(__file__))+"/flair/flair.py", help='flair path')
parser.add_argument('-g', '--genome', action='store', dest='g', default=os.path.dirname(os.path.realpath(__file__))+"/GRCh38.primary_assembly.genome.fa", help='path to genome')
#parser.add_argument('-x', '--minimap', action='store', dest='x', default="/private/groups/brookslab/bin/minimap2", help='path to minimap')
parser.add_argument('-k', '--remapSize', action='store', dest='k', default=0, type=int, help='size of area to remap - only remaps if this is specified')
parser.add_argument('-t', '--transcriptome', action='store', dest='t', default=os.path.dirname(os.path.realpath(__file__))+"/gencode.v37.annotation.gtf", help='path to transcriptome (.gtf)')
parser.add_argument('-n', '--spliceJunctions', action='store', dest='n', default=os.path.dirname(os.path.realpath(__file__))+"/intropolis.liftover.hg38.junctions.sorted.txt", help='path to splice junction file (.txt)')
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
	tempsamflags = {}
	if args.s:
		print('making sam file')
		process = subprocess.Popen('samtools view -h -o ' + correctQ + '.sam ' + correctQ + '.bam',stdout=subprocess.PIPE, shell=True)
		print(process.communicate()[0].strip())
	correctQ = prefix if os.path.exists(prefix + '.sam') else prefix + '.aligned'
	# allMatches = set(allMatches)
	readLength = {}
	print('processing sam file')
	sam = open(correctQ + '.sam', 'r')
	flagbinary = {'0':0, '2048':0, '16':1, '2064':1, '256':0, '272':1}
	mappingLocs = {}
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
				# print(gene, loc, locchr, line[:4])
				if locchr == line[2] or currGene == None:
					if abs(int(line[3]) - loc) < currDist:
						currDist, currGene = abs(int(line[3]) - loc), gene
			index = new_fusions_found[readToFusion[line[0]]]['readNames'].index(line[0])
			if flagbinary[line[1]] == 1:#!= flagbinary[mappingLocs[readToFusion[line[0]]][line[0]][0][-1]]:
				temp = locs[0]
				locs[0] = readLength[line[0]] - locs[1]
				locs[1] = readLength[line[0]] - temp
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[0], new_fusions_found[readToFusion[line[0]]][currGene]['right'][index], currGene])
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[1], new_fusions_found[readToFusion[line[0]]][currGene]['left'][index], currGene])
			else:
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[0], new_fusions_found[readToFusion[line[0]]][currGene]['left'][index], currGene])
				mappingLocs[readToFusion[line[0]]][line[0]].append([locs[1], new_fusions_found[readToFusion[line[0]]][currGene]['right'][index], currGene])

###TODO: fastqdist, ssdist, fastqcov, genecov, tooclose, ismito, readsup
	print('filtering fusions', len(mappingLocs.keys()))
	filteredFusions = []
	for f in mappingLocs:
		if len(mappingLocs[f].keys()) >= args.l:
			# print(f, mappingLocs[f]
			c = 0
			fastqcov, genecov, tooshort, genecovfinal, shortmap = [], {}, False, [0], 0
			for l in list(f):
				new_fusions_found[f][l]['chr'] = max(set(new_fusions_found[f][l]['chr']), key=new_fusions_found[f][l]['chr'].count)
				genecov[l] = []
				for i in range(len(new_fusions_found[f][l]['left'])):
					genecov[l].append([new_fusions_found[f][l]['right'][i]-new_fusions_found[f][l]['left'][i], new_fusions_found[f][l]['left'][i], new_fusions_found[f][l]['right'][i]])
			for l in genecov:
				genecov[l].sort()
				# if median([x[0] for x in genecov[l][:int(len(genecov[l])/4)]]) < 200: tooshort = True
				midcov, bestcov = genecov[l][int(len(genecov[l])/2)], 0
				if midcov[0] < 100: tooshort, shortmap = True, midcov[0]
				if not (l[:3] == 'chr' and '-' in l) and l in allGeneLoc:
					for t in allGeneLoc[l][1:]:
						if midcov[1] >= t[1]-50 and midcov[2] <= t[2]+50:
							if float(midcov[0])/(t[2]-t[1]) > bestcov: bestcov = float(midcov[0])/(t[2]-t[1])
					genecovfinal.append(bestcov)
			if max(genecovfinal) > 0.95:
				metadata.append(['--'.join(list(f)), 'geneCov', max(genecovfinal), list(mappingLocs[f].keys())])
				continue
			if tooshort:
				metadata.append(['--'.join(list(f)), 'tooShortMapping', shortmap, list(mappingLocs[f].keys())])
				continue
			maxlen = 0
			for r in mappingLocs[f]:
				if len(mappingLocs[f][r]) > maxlen: maxlen = len(mappingLocs[f][r])
			possPromLocs = [[] for i in range(maxlen)]
			tempLocs = [[] for i in range(maxlen)]
			firstgene = max(set([mappingLocs[f][x][0][2] for x in mappingLocs[f]]), key=[mappingLocs[f][x][0][2] for x in mappingLocs[f]].count)
			for r in mappingLocs[f]:
				print(f, r, mappingLocs[f][r])
				tempfastqcov = 0
				for l in range(0, len(mappingLocs[f][r]), 2):
					tempfastqcov += mappingLocs[f][r][l+1][0]-mappingLocs[f][r][l][0]
				fastqcov.append(float(tempfastqcov)/readLength[r])
				mappingLocs[f][r].sort()
				# if f == frozenset(['SLC25A24', 'NBPF6']): print(r)
				# if f == frozenset(['SLC25A24', 'NBPF6']): print(mappingLocs[f][r])
				if mappingLocs[f][r][0][2] != firstgene:
					mappingLocs[f][r] = mappingLocs[f][r][::-1]
					for x in range(len(mappingLocs[f][r])): mappingLocs[f][r][x][0] = readLength[r] - mappingLocs[f][r][x][0]
				for x in range(1, len(mappingLocs[f][r])-2, 2):
					if mappingLocs[f][r][x][2] != mappingLocs[f][r][x-1][2]:
						temp = mappingLocs[f][r][x]
						mappingLocs[f][r][x] = mappingLocs[f][r][x+1]
						mappingLocs[f][r][x+1] = temp
				# if f == frozenset(['SLC25A24', 'NBPF6']): print(mappingLocs[f][r])
				myrange = list(range(int(len(mappingLocs[f][r])/2))) + list(range(-int(len(mappingLocs[f][r])/2), 0)) if len(mappingLocs[f][r]) < maxlen else range(len(mappingLocs[f][r]))
				for l in myrange:
					possPromLocs[l].append(mappingLocs[f][r][l])
					# tempLocs[l].append(mappingLocs[f][r][l] + [r, new_fusions_found[f][mappingLocs[f][r][l][2]][mappingLocs[f][r][l][1]][new_fusions_found[f]['readNames'].index(r)]])
			if median(fastqcov) < 0.8:
				metadata.append(['--'.join(list(f)), 'fastqCov', median(fastqcov), list(mappingLocs[f].keys())])
				continue
			distToProm, gene5, geneorder = 1000000000000000000000, None, []
			fastqdist, ssdist, tooclose, closedist = [], [], False, 0
			isMito, mitogene= False, None



if int(args.k) > 0:
	readNames = []
	print('remapping reads to fusions')
	with open(outfilename + "readNames.txt", 'r') as names:
		for line in names:
			readNames.append(line.rstrip())
	readNames = set(readNames)
	#Filter reads to only the double mapped reads
	leadingChar = '@' if args.r.split('.')[-1] == 'fastq' or args.r.split('.')[-1] == 'fq' else '>'
	with open(args.r, 'r') as reads, open(outfilename + "Filtered.fa", "w") as faOut:
		writeRead = False
		c, d = 0, 0
		for line in reads:
			c += 1
			if c % 20000000 == 0: print(c, d)
			if line[0] == leadingChar:
				if line.rstrip('\n').lstrip(leadingChar) in readNames:
					writeRead = True
					d += 1
					if leadingChar == '@': faOut.write(">" + line.lstrip('@'))
					else: faOut.write(line)
				else:
					writeRead = False
			elif writeRead:
				faOut.write(line)
				writeRead = False
	print('reads filtered')
	fusions = {}
	firstLine = []
	with open(outfilename + "Fusions.tsv", 'r') as thesefusions, open(outfilename + 'Locs.bed', 'w') as bedFile:
		for line in thesefusions:
			line = line.rstrip().split('\t')
			if line[0][0] != '#' and len(line) > 6:
				chr1, center1 = line[5].split('-')[-3:-1]
				chr2, center2 = line[6].split('-')[-3:-1]
				if int(center1) > args.k and int(center2) > args.k:
					name1 = line[0] + '->' + '-'.join(line[5].lstrip("3'-").split('-')[:-2])
					name2 = line[0] + '->' + '-'.join(line[6].lstrip("5'-").split('-')[:-2])
					fusions[line[0]] = {'line':line}
					fusions[line[0]]['-'.join(line[5].lstrip("3'-").split('-')[:-2])] = {'side':"3'", 'chr':chr1, 'bp':int(center1), 'left':[], 'right':[], 'reads':[], 'mapQ':[]}
					fusions[line[0]]['-'.join(line[6].lstrip("5'-").split('-')[:-2])] = {'side':"5'", 'chr':chr2, 'bp':int(center2), 'left':[], 'right':[], 'reads':[], 'mapQ':[]}
					bedFile.write('\t'.join([chr1, str(int(center1)-args.k), str(int(center1) + args.k), name1]) + '\n')
					bedFile.write('\t'.join([chr2, str(int(center2)-args.k), str(int(center2) + args.k), name2]) + '\n')
			elif line[0][0] == '#':
				firstLine = line
	process = subprocess.Popen('bedtools getfasta -fi ' + args.g + ' -bed ' + outfilename + 'Locs.bed' + ' -fo ' + outfilename + 'Genome.fa' + ' -name; ' +
							   'minimap2 -a ' + outfilename + 'Genome.fa ' + outfilename + "Filtered.fa" + ' > ' + outfilename + 'Remapped.sam; ' +
							   "sam2bed < " + outfilename + 'Remapped.sam' + ' > ' + outfilename + 'Remapped-unfilt.bed',stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())
	maxMapQ = 0
	with open(outfilename + 'Remapped-unfilt.bed', 'r') as remapped:
		for line in remapped:
			line = line.rstrip().split('\t')
			if int(line[4]) > maxMapQ: maxMapQ = int(line[4])
			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['mapQ'].append(int(line[4]))
			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['reads'].append(line[3])
			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['left'].append(int(line[1])-args.s)
			fusions[line[0].split('->')[0]][line[0].split('->')[-1]]['right'].append(int(line[2])-args.s)
	newFusions = open(outfilename + 'FusionsRemapped.tsv', 'w')
	firstLine.insert(2, 'confirmed reads')
	newFusions.write('\t'.join(firstLine) + '\n')
	doubleMappedReads = []
	for fusion in fusions:
		good = 0
		fusionReads = []
		mapScores = []
		for loc in fusions[fusion]:
			if loc not in ['line', 'conf reads']:
				if len(fusions[fusion][loc]['left']) > 0:
					leftAvg, rightAvg = sum(fusions[fusion][loc]['left'])/len(fusions[fusion][loc]['left']), \
										sum(fusions[fusion][loc]['right'])/len(fusions[fusion][loc]['right'])
					if abs(leftAvg + rightAvg) > 50:
						good += 1
						# if abs(leftAvg) > abs(rightAvg):
						# 	fusions[fusion][loc]['bp'] += max(fusions[fusion][loc]['right'])
						# else:
						# 	fusions[fusion][loc]['bp'] += min(fusions[fusion][loc]['left'])
						fusionReads.append(fusions[fusion][loc]['reads'])
						mapScores += fusions[fusion][loc]['mapQ']
						# if fusions[fusion][loc]['side'] == "3'":
						# 	fusions[fusion]['line'][5] = '-'.join(fusions[fusion]['line'][5].split('-')[:-1] + [str(fusions[fusion][loc]['bp'])])
						# else:
						# 	fusions[fusion]['line'][6] = '-'.join(fusions[fusion]['line'][6].split('-')[:-1] + [str(fusions[fusion][loc]['bp'])])
		if good >= 2:
			shared = list(set(fusionReads[0]) & set(fusionReads[1]))
			if len(shared) > 0:#good >= 2:
				for i in shared: doubleMappedReads.append(fusion + i)
				fusions[fusion]['line'].insert(2, str(len(shared)))
				fusions[fusion]['line'][3] = str(round((sum(mapScores)/len(mapScores))/float(maxMapQ), 3))
				newFusions.write('\t'.join(fusions[fusion]['line']) + '\n')
	newFusions.close()
	doubleMappedReads = set(doubleMappedReads)
	goodReadNames = []
	with open(outfilename + 'Remapped-unfilt.bed', 'r') as remapped, open(outfilename + 'Remapped.bed', 'w') as remapFilt:
		for line in remapped:
			line = line.split('\t')
			if line[0].split('->')[0] + line[3] in doubleMappedReads:
				remapFilt.write('\t'.join(line))
				goodReadNames.append(line[3])
	readsOut = open(outfilename + 'RemappedReads.bed', 'w')
	goodReadNames = set(goodReadNames)
	for line in open(outfilename + "Reads.bed", 'r'):
		temp = line.split('\t')
		if temp[3].split('-.-')[1] in goodReadNames:
			readsOut.write(line)


	process = subprocess.Popen('python3 ' + os.path.dirname(os.path.realpath(__file__)) + '/makeAlnSeq.py -f ' + outfilename +
							   'Genome.fa -r ' + outfilename + 'Remapped-unfilt.bed; rm ' + outfilename + 'Remapped-unfilt.bed',stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())


if not args.i:
	left, right = [], []
	for line in open(outfilename + 'Fusions.tsv', 'r'):
		if line[0] != '#':
			line = line.rstrip().split('\t')
			fusion = sorted(line[0].split('--'))
			# print(fusion)
			left.append(fusion[0])
			right.append(fusion[1])
	left = set(left)
	right = set(right)
	# print("left")
	# print(left)
	leftOut = open(outfilename + 'Reads-l.bed', 'w')
	rightOut = open(outfilename + 'Reads-r.bed', 'w')
	leftReads, rightReads = [], []
	c = 0
	for line in open(outfilename + 'Reads.bed', 'r'):
		line2 = line.strip().split('\t')
		geneName = line2[3].split('-.-')[2].split('/')[0]
		if geneName in left:
			c += 1
			leftOut.write(line)
		elif geneName in right:
			rightOut.write(line)
	leftOut.close()
	rightOut.close()
	# print(c)
	#print('python3 ' + args.f + ' collapse -g ' + args.g + ' -r ' + args.r + ' --generate_map -q ' + outfilename + 'Reads.bed -o ' + prefix + '.fusions',)
	#print('python3 ' + args.f + ' collapse -g ' + args.g + ' -r ' + args.r + ' --generate_map -q ' + outfilename + 'Reads-l.bed -o ' + prefix + '.fusions.l')
	process = subprocess.Popen(
		#collapse breaks.simplen/he_v2.3.5_pass.fastq --generate_map -q 17-03-2021he_v2.3.5_passReads-1.bed -o he.fusions.collapse
		'python3 ' + args.f + ' collapse --stringent --temp_dir temp_dir_l -g ' + args.g + ' -r ' + args.r + ' --generate_map -q ' + outfilename + 'Reads-l.bed -o ' + prefix + '.fusions.l' +
		'; python3 ' + args.f + ' collapse --stringent --temp_dir temp_dir_r -g ' + args.g + ' -r ' + args.r + ' --generate_map -q ' + outfilename + 'Reads-r.bed -o ' + prefix + '.fusions.r',
		stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())

if not args.j and os.path.exists(prefix + '.fusions.l.isoform.read.map.txt'):
	fusionisos = {}
	fusionHeader = ''
	for line in open(outfilename + 'Fusions.tsv', 'r'):
		if line[0] == '#': fusionHeader = line
		else:
			line = line.strip().split('\t')
			temp = line[0].split('--')
			fusionisos[line[0]] = {'reads':[], 'isos':[], 'line':line}
	for line in open(outfilename + 'Reads.bed', 'r'):
		line = line.split('\t')
		info = line[3].split('-.-')
		if info[0] in fusionisos.keys():
			fusionisos[info[0]]['reads'].append(info[1])
	lIsoSup, rIsoSup = {}, {}
	lIsoKey, rIsoKey = {}, {}
	c = 0
	for line in open(prefix + '.fusions.l.isoform.read.map.txt', 'r'):
		temp = line.strip().split('\t')
		lIsoKey['lIso' + str(c)] = temp[0]
		c += 1
		for i in temp[1].split(','):
			lIsoSup[i] = temp[0]#'lIso' + str(c) #temp[0]
	c = 0
	for line in open(prefix + '.fusions.r.isoform.read.map.txt', 'r'):
		temp = line.strip().split('\t')
		rIsoKey['rIso' + str(c)] = temp[0]
		c += 1
		for i in temp[1].split(','):
			rIsoSup[i] = temp[0]#'rIso' + str(c) #temp[0]
	lSet, rSet = set(lIsoSup.keys()), set(rIsoSup.keys())
	for i in fusionisos.keys():
		fusionisos[i]['reads'] = set(fusionisos[i]['reads'])
		for j in fusionisos[i]['reads']:
			inl, inr = '.', '.'
			# TESTING
			if j in lSet:
				if lIsoSup[j].split('-.-')[0]==i:
					inl = lIsoSup[j]
			if j in rSet:
				if rIsoSup[j].split('-.-')[0]==i:
					inr = rIsoSup[j]
			if inl != '.' or inr != '.':
				fusionisos[i]['isos'].append('-=-'.join([inl, inr]))
			# END TESTING
	c = 0
	#print(lSet)
	isoToCombName = {}
	fusionsOut = open(outfilename + 'IsoformFusions.tsv', 'w')
	fusionsOut.write(fusionHeader)
	#print(fusionisos)
	for i in fusionisos.keys():
		combIsoCounts = Counter(fusionisos[i]['isos'])
		#if c < 30: print(i, combIsoCounts)
		totIsos, doubleIsos, readsSupIsos = 0, 0, 0
		for j in combIsoCounts.keys():
			temp = j.split('-=-')
			numsides = 'single'
			if temp[0] != '.' and temp[1] != '.':
				doubleIsos += 1
				numsides = 'double'
			totIsos += 1
			#print(tempr)
			readsSupIsos += combIsoCounts[j]
			if temp[0] != '.':
				if temp[0].split('-.-')[1] not in isoToCombName: isoToCombName[temp[0].split('-.-')[1]] = []
				isoToCombName[temp[0].split('-.-')[1]].append("#" + i)
				isoToCombName[temp[0].split('-.-')[1]].append(str(combIsoCounts[j]) + '|combIso' + str(c) + '|' + numsides)
			if temp[1] != '.':
				if temp[1].split('-.-')[1] not in isoToCombName: isoToCombName[temp[1].split('-.-')[1]] = []
				isoToCombName[temp[1].split('-.-')[1]].append("#" + i)
				isoToCombName[temp[1].split('-.-')[1]].append(str(combIsoCounts[j]) + '|combIso' + str(c) + '|' + numsides)
			c += 1
		fusionisos[i]['line'][1] = str(doubleIsos) + '/' + str(totIsos) + '/' + str(readsSupIsos) + '/' + fusionisos[i]['line'][1]
		if totIsos > 0:
			fusionsOut.write('\t'.join(fusionisos[i]['line']) + '\n')
	readsOut = open(outfilename + 'IsoformReads.bed', 'w')
	#TESTING
	# multiGeneIsos = []
	# for i in isoToCombName.keys():
	# 	fusionList = []
	# 	for j in isoToCombName[i]:
	# 		if j[0] == '#': fusionList.append(j)
	# 	if len(Counter(fusionList).keys()) > 1:
	# 		multiGeneIsos.append([i, fusionList])
	# print(multiGeneIsos)
	# print(len(multiGeneIsos))
	c = 0
	written = []
	for i in ['l', 'r']:
		for line in open(prefix + '.fusions.' + i + '.isoforms.bed', 'r'):
			line = line.strip().split('\t')
			info = line[3].split('-.-')
			if info[1] in isoToCombName.keys():
				# c += 1
				# if c < 10:
				# 	print(isoToCombName[info[1]])
				for j in isoToCombName[info[1]]:
					if j[0] != '#':
						temp = info
						info[1] = j
						line[3] = '-.-'.join(info)
						if line[3] not in written:
							written.append(line[3])
							readsOut.write('\t'.join(line) + '\n')
	fusionsOut.close()
	readsOut.close()
	process = subprocess.Popen('rm -f ' + prefix + '.fusions.l.isoforms.fa ' + prefix + '.fusions.l.isoforms.bed ' + prefix + '.fusions.l.isoform.read.map.txt' + prefix + '.fusions.r.isoforms.fa ' + prefix + '.fusions.r.isoforms.bed ' + prefix + '.fusions.r.isoform.read.map.txt ' + outfilename + 'Reads-l.bed ' + outfilename + 'Reads-r.bed',stdout=subprocess.PIPE, shell=True)
	print(process.communicate()[0].strip())


# read on one side can be matched to multiple reads on the other (in combIsoCounts)
# must pick best
# why are reads from diff. fusions getting intp fusion?
# make isoform support stats lign up
# when multiple fusions in the same gene, it messes up
# check for any fusions with 2 unique loci - see if they look the same on the og version of testing section and new version