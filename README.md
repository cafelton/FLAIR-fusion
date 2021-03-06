# FLAIR-fusion

Requires FLAIR (https://github.com/BrooksLabUCSC/flair), python3 and numpy

minimap2 and bedtools2.28 must be in your path

Also download intropolis.liftover.hg38.junctions.sorted.txt at https://drive.google.com/file/d/10Kz7lzVQlNF2ANoEKLcYIXPgfRubxQCQ/view?usp=sharing
and save it to the same folder as 19-03-2021-fasta-to-fusions-pipe.py


First: either download gencode.v37.annotation-short.gtf if using GRCh38 from https://drive.google.com/file/d/1oEUrrom8evGk9b1m7CSCp0PKdlEbWq9u/view?usp=sharing

or run makeShortAnno with

python makeShortAnno.py /other-folder/gene-annotation.gtf

Next: run the full pipeline with

python3 19-03-2021-fasta-to-fusions-pipe.py -r file.fastq -f path/to/flair.py -g /path/to/genome.fa -t /path/to/anno.gtf -a /path/to/anno-short.gtf


Required (run in python3)

-r --reads    fastq of fasta file of long reads (nanopore or pacbio)

-f --flair    path to flair

-g --genome    path to genome .fa

-t --transcriptome path to gene annotation .gtf

-a --anno    short gene annotation file.gtf from makeShortAnno


Optional


-o --output   output prefix (added to fastq prefix) default-today's date

-b --buffer length of buffer for combining nearby regions and determining distinct loci. default 50000


-s --samConvert    whether to take .bam and convert to .sam (True = convert .bam (from fq prefix) to .sam) - not necessary if you're doing the alignment step

-y --includeMito    whether to include fusions that are in the mitochondria (True=include) default=False not reccommended


-k --remapSize size of area around breakpoint to remap  default-0, reccommended-500

-i --callIsoforms    whether to detect fusion isoforms (True=already detected or don't want to detect, dont run) default=False

-j --matchFusionIsos    whether to match isoforms to fusions (True=already matched or dont want to match, dont run) default=False


-d --detectFusions   whether to detect fusions (True=already detected, dont run) default=False


-p --bedProcess whether to align and correct reads (True=I already have a processed .bed file with the filename in the form fastqName-bedtools-genes-short.bed). default = False

-u --flairAlign    whether to align reads (True=I already have an aligned .bed file (name stored in -m))  default=False

-c --flairCorrect    whether to correct reads (True=I already have a corrected .bed file (name stored in -m)  default=False

-m --bedfile    name of aligned.bed file or corrected.bed file if -u or -c is selected


Output:

if detecting fusions: one args.o-fastq-Fusions.tsv file and one args.o-fastq-Reads.bed file with only chimeric reads

if remapping (-k > 0): extra args.o-fastq-FusionsRemapped.tsv and args.o-fastq-Remapp-seq.bed file. The .bed file will not be in standard form, since we make synthetic chromosomes around the fusion breakpoints.

if detecting isoforms: extra args.o-fastq-IsoformFusions.tsv and args.o-fastq-IsoformsReads.bed file with chimeric isoforms. These isoforms will be renamed and not the names of any of your reads.
