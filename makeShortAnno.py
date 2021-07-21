import sys

#run with python makeShortAnno.py annotation.gtf
out = open('.'.join(sys.argv[1].split('.')[:-1]).split('/')[-1] + "-short.gtf", 'w')
for line in open(sys.argv[1]):
    if line[0] == '#':
        out.write(line)
    else:
        temp = line.split('\t')
        if temp[2] == 'gene':
            out.write(line)

