.PHONY: all
all: src/htslib-1.21 data/ce11.fa data/hg38_long_mrna.fa
src/htslib-1.21:
	env -C src wget -4 \
		https://sourceforge.net/projects/samtools/files/samtools/1.21/htslib-1.21.tar.bz2/download \
		-O htslib-1.21.tar.bz2
	env -C src tar xvjf htslib-1.21.tar.bz2

data/ce11.fa:
	curl https://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz | gunzip > data/ce11.fa

data/hg38_long_mrna.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz | gunzip | seqtk seq -L 500 > data/hg38_long_mrna.fa

