import random

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(seq):
    bases = list(seq)
    # bases = reversed([dcomplement.get(base,base) for base in bases])
    bases = [complement[x] for x in bases]
    bases = bases[::-1]
    bases = ''.join(bases)
    return bases

transcripts = {}

# print(genome['SIRV1'])
#dont'f forget to add directionality - reverse complement function?
genes = {}
##format: genes[geneName] = {chr:'chrZ', 'start':348756, 'end':28457246,
        # 'isoformName1': ['+/-', [exon1start, exon1end], [exon2start, exon2end]],
        #'isoformName2': ['+/-', [exon1start, exon1end], [exon2start, exon2end]]}
##Example: genes[SGIP1]: = {chr:'chr1', 'start':66531592, 'end':66753139,
        # 'SGIP1-213' ['+', [66533604, 66533672], [66534245, 66534368], [66590469, 66590606], [66625846, 66626043]]
        # 'SGIP1-205 ['+', [66533614, 66533672], [66534245, 66534368], [66625846, 66625910], [66633069, 66633094], [66635943, 66636015], [66639776, 66639833], [66642809, 66642864], [66643543, 66643719], [66670994, 66671019], [66671943, 66671995], [66673280, 66673366], [66677003, 66677096], [66679677, 66679742]]}

for line in open("/private/groups/brookslab/reference_annotations/gencode.v37.annotation.gtf"):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'exon' and line[8].split('; gene_type "')[1].split('"')[0] == "protein_coding":
            gene = line[8].split('; gene_name "')[1].split('"')[0]
            name = line[8].split('; transcript_name "')[1].split('"')[0]
            if 'orf' not in gene:
                if gene not in genes.keys():
                    genes[gene] = {'chr':line[0], 'start':int(line[4]), 'end':int(line[3]), 'dir':line[6]}
                    # leftG[gene] = []
                    # rightG[gene] = []
                if name not in genes[gene].keys():
                    genes[gene][name] = [line[6]]
                if line[6] == '+':
                    genes[gene][name].append([int(line[3]), int(line[4])])
                    if genes[gene]['end'] < int(line[3]): genes[gene]['end'] = int(line[3])
                else:
                    genes[gene][name] = [genes[gene][name][0], [int(line[3]), int(line[4])]] + genes[gene][name][1:]
                    if genes[gene]['start'] > int(line[4]): genes[gene]['start'] = int(line[4])
#print(transcripts)
print('gtf processed')


###Overview of simulation method from this paper: https://www.sciencedirect.com/science/article/pii/S2405471221002076?via%3Dihub#sec6
# Fusion transcripts were simulated by randomly selecting two protein coding transcripts (using Ensembl v75 annotation).
# For each of the selected transcripts, a random number of exons are used to generate the fusion. If the transcript was
# selected as the donor, then the number of exons incorporated are counted from the beginning of the transcript.
# Conversely, if the transcript was selected as the acceptor, then the number of exons incorporated are counted from
# the end of the transcript. Selected transcripts are fused only at the exon-intron boundaries. Using the exon coordinates
# for each selected transcript, a synthetic fusion sequence is generated using the GRCh37.75 genome and biopython.
# A reference index is generated for the synthetic fusion sequence using RSEM v1.2.31 with the STAR 2.4.2a aligner.


out = open("gencode37-fusion-sim-test-06-08-2023.fasta", 'w')
out2 = open('gencode37-fusion-sim-test-06-08-2023-truth.tsv', 'w')
out3 = open('gencode37-fusion-sim-test-06-08-2023-truth.bed', 'w')
count = 0
leftG, rightG = {}, {}
isos = {}
### format: leftG[gene] = ['genomic sequence(ACGT) of all exons up to cutpoint in isoform1', 'genomic sequence(ACGT) of all exons up to cutpoint in isoform2']
##Example: leftG[TEX26] =  ['ATACTTGAATTGTTCAATTCTGAAAGAATATTGGAACCTCCCACAATAAATATATTTTTATCAA', 'GAAGAAAATCTGAAAATCAATGGAAGCACGAGGACATTCCTAGTCATTTTCTCAATTATCAAGGAAAAATAAGATGCAAATAGCTTCAA']
while count < 100:
    g = random.choice(list(genes.keys()))
    if len(genes[g]) > 5 and genes[g]['end']>genes[g]['start']:
        genes[g]['cutpoint'] = random.randrange(genes[g]['start'], genes[g]['end'])
        written = False
        for t in genes[g]:
            if t not in ['chr', 'start', 'end', 'cutpoint', 'dir']:
                # print(g, genes[g])
                # print(genes[g][t])
                if genes[g][t][1][1] < genes[g]['cutpoint'] and genes[g][t][-1][0] > genes[g]['cutpoint']: #check if cutpoint is after first exon
                    leftT, rightT = [], []#"", ""
                    leftCut, rightCut = 0, 0
                    for i in genes[g][t][1:]:
                        if i[1] < genes[g]['cutpoint']:
                            leftT.append((i[0],i[1])) #+= genome[genes[g]['chr']][i[0]:i[1]]
                            leftCut = i[1]
                        elif i[0] >= genes[g]['cutpoint']:
                            rightT.append((i[0], i[1])) #+= genome[genes[g]['chr']][i[0]:i[1]]
                            if rightCut == 0: rightCut = i[0]
                    if g not in leftG: leftG[g] = []
                    if g not in rightG: rightG[g] = []
                    if g not in isos: isos[g] = []
                    leftG[g].append(leftT)
                    rightG[g].append(rightT)
                    isos[g].append([t, leftCut, rightCut])
                    written = True
                    # if count < 5: print(g, t, genes[g][t])
        if written:
            count += 1


genome = {}
last = None
for line in open("/private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa"):
    if line[0] == '>':
        if last != None:
            genome[last] = ''.join(genome[last])
             # print(genome[last][-300:])
        last = line.rstrip('\n').split()[0].lstrip('>')
        # print(last)
        genome[last] = []
    else:
        genome[last].append(line.rstrip('\n'))
print('genome processed')






done = []
count = 0
out2.write('gene1\tgene2\tnumOfFusionIsos\tgene1chr\tgene1breakpoint\tgene1isos\tgene2chr\tgene2breakpoint\tgene2isos\n')
for g1 in leftG:
    for g2 in leftG:
        if g1 != g2 and g1 not in done and g2 not in done:
            c = 0
            if genes[g1]['dir'] == genes[g2]['dir']:
                # if count < 3: print(g1, leftG[g1], g2, rightG[g2])
                minIsos = min([len(leftG[g1]), len(rightG[g2])])
                if genes[g1]['dir'] == '+':
                    out2.write('\t'.join([g1, g2,
                                          str(minIsos), genes[g1]['chr'], genes[g1]['dir'], str(genes[g1]['cutpoint']),
                                          ','.join([x[0] + '#' + str(x[1]) for x in isos[g1][:minIsos]]), genes[g2]['chr'], genes[g2]['dir'], str(genes[g2]['cutpoint']),
                                          ','.join([x[0] + '#' + str(x[2]) for x in isos[g2][:minIsos]])]) + '\n')
                else:
                    out2.write('\t'.join([g2, g1, str(minIsos), genes[g2]['chr'], genes[g2]['dir'], str(genes[g2]['cutpoint']),
                                          ','.join([x[0] + '#' + str(x[2]) for x in isos[g2][:minIsos]]), genes[g1]['chr'], genes[g1]['dir'], str(genes[g1]['cutpoint']),
                                          ','.join([x[0] + '#' + str(x[1]) for x in isos[g1][:minIsos]])]) + '\n')
                # print(g1 + '\t' + g2 + '\t' + str(minIsos))
                for i in range(minIsos):
                    out.write('>jointest' + '-' + isos[g1][i][0] + '-' + isos[g2][i][0] + '-' + str(count) + '-' + str(c) + '\n')
                    leftSeq, rightSeq = '',''
                    for exon in leftG[g1][i]:
                        leftSeq += genome[genes[g1]['chr']][exon[0]:exon[1]]
                    for exon in rightG[g2][i]:
                        rightSeq += genome[genes[g2]['chr']][exon[0]:exon[1]]
                    if genes[g1]['dir'] == '+':
                        out.write(leftSeq + rightSeq + '\n') #leftG[g1][i] + rightG[g2][i] + '\n')
                    else:
                        out.write(reverse_complement(leftSeq + rightSeq) + '\n') #leftG[g1][i] + rightG[g2][i]) + '\n')
                    c += 1
                #write to truth bed file
                for i in range(minIsos):
                    isoName = isos[g1][i][0] + '--' + isos[g2][i][0] if genes[g1]['dir'] == '+' else isos[g2][i][0] + '--' + isos[g1][i][0]
                    leftBed = [genes[g1]['chr'], leftG[g1][i][0][0], isos[g1][i][1], isoName, 60, genes[g1]['dir'], leftG[g1][i][0][0], isos[g1][i][1], '255,0,0', len(leftG[g1][i]), [], []]
                    for exon in leftG[g1][i]:
                        leftBed[-2].append(exon[1]-exon[0])
                        leftBed[-1].append(exon[0]-leftG[g1][i][0][0])
                    leftBed[-2] = ','.join([str(x) for x in leftBed[-2]])
                    leftBed[-1] = ','.join([str(x) for x in leftBed[-1]])

                    rightBed = [genes[g2]['chr'], isos[g2][i][2], rightG[g2][i][-1][1],
                               isoName, 60, genes[g2]['dir'], isos[g2][i][2], rightG[g2][i][-1][1], '255,0,0', len(rightG[g2][i]), [], []]
                    for exon in rightG[g2][i]:
                        rightBed[-2].append(exon[1] - exon[0])
                        rightBed[-1].append(exon[0] - isos[g2][i][2])
                    rightBed[-2] = ','.join([str(x) for x in rightBed[-2]])
                    rightBed[-1] = ','.join([str(x) for x in rightBed[-1]])

                    out3.write('\t'.join([str(x) for x in leftBed]) + '\n')
                    out3.write('\t'.join([str(x) for x in rightBed]) + '\n')

                count += 1

            # c = 0
            # for i in range(min([len(leftG[g2]), len(rightG[g1])])):
            #     out.write('>jointest' + '-' + g1 + '-' + g2 + '-' + str(count) + '-' + str(c) + '\n')
            #     out.write(leftG[g2][i] + rightG[g1][i] + '\n')
            #     c += 1
            # count += 1
                done.append(g1)
                done.append(g2)
print(count)
print('fusions done')
# count = 0
# for x in range(8000):
#     g = random.choice(list(genes.keys()))
#     c = 0
#     while c < random.randint(0,10) or c < len(genes[g].keys())-4:
#         t = random.choice(list(genes[g].keys()))
#         if t not in ['chr', 'start', 'end', 'cutpoint']:
#             out.write('>test' + str(count) + '\n')
#             thisT = ""
#             for i in genes[g][t][1:]:
#                 thisT += genome[genes[g]['chr']][i[0]:i[1]]
#             out.write(thisT + '\n')
#             count += 1
#             c += 1
# out.close()
# print('other genes done')