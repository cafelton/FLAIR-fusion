import argparse
parser = argparse.ArgumentParser(description='fusion caller parse options', usage='python id-fusions.py -s genomeAlign.sam -r reads.fa -o outputPrefix OR -g genomeAlign.psl -t transcriptAlign.psl -o outputPrefix')
parser.add_argument('-i', '--reads', action='store', dest='i', default="", help='input.txt')
args = parser.parse_args()

prev = [""]*12
with open(args.i, 'r') as inFile, open('.'.join(args.i.split('.')[:-1]) + '-short.bed', "w") as out:
	for line in inFile:
		line = line.strip().split('\t')
		if len(line) > 20:
			if line[:6] != prev[:6]:
				if line[20] != '.':
					name = line[20]#.strip().split('; gene_name "')[1].split('"')[0]
					#line[3] = line[3] + '--' + name
					line.append(name)
				else:
					#line[3] = line[3] + '--' + line[0] + '-' + str(round(int(line[1]), -4))
					line.append(line[0] + '-' + str(round(int(line[1]), -4)))
				prev[3] = prev[3] + '--' + prev[-1].split('...')[0]
				if len(prev[0]) <= 3:
					prev[0] = 'chr' + prev[0]
				out.write("\t".join(prev[:12]) + '\n') #line[:11]
				prev = line
			if line[14] == 'gene':
				name = line[20].strip().split('; gene_name "')[1].split('"')[0]
				if '/' in prev[-1]:
					temp = round(float(line[21])/(int(line[16])-int(line[15])), 5)
					if (line[20].split('; gene_type "')[1].split('"')[0] == 'protein_coding' and prev[-1].split('...')[2] == 'protein_coding') or \
						(line[20].split('; gene_type "')[1].split('"')[0] != 'protein_coding' and prev[-1].split('...')[2] != 'protein_coding'):
						if temp > float(prev[-1].split('...')[1]) and temp < 1 and (int(line[16])-int(line[15])) > (int(line[2]) - int(line[1])):
							prev[-1] = "/".join([name, line[18], line[15], line[16]]) + '...' + str(temp) + '...' + line[20].split('; gene_type "')[1].split('"')[0]
					elif line[20].split('; gene_type "')[1].split('"')[0] == 'protein_coding':
						prev[-1] = "/".join([name, line[18], line[15], line[16]]) + '...' + str(temp) + '...' + line[20].split('; gene_type "')[1].split('"')[0]
				else:
					prev[-1] = "/".join([name, line[18], line[15], line[16]]) + '...' + str(round(float(line[21])/(int(line[16])-int(line[15])), 5)) + \
							   '...' + line[20].split('; gene_type "')[1].split('"')[0]
