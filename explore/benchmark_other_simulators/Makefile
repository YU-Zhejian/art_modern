.PHONY: all
all: src/htslib-1.21 src/gsl-2.8 \
	data/ce11.fa data/hg38_long_mrna.fa \
	data/yeast.fa data/e_coli.fa data/soybean.fa data/lambda_phage.fa

src/htslib-1.21:
	env -C src wget -4 \
		https://sourceforge.net/projects/samtools/files/samtools/1.21/htslib-1.21.tar.bz2/download \
		-O htslib-1.21.tar.bz2
	env -C src tar xvjf htslib-1.21.tar.bz2

src/gsl-2.8:
	env -C src wget -4 \
    		https://ftpmirror.gnu.org/gsl/gsl-2.8.tar.gz \
    		-O gsl-2.8.tar.gz
	env -C src tar xzvf gsl-2.8.tar.gz

data/ce11.fa:
	curl https://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz | gunzip > data/ce11.fa

data/hg38_long_mrna.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz | \
	gunzip | \
	seqtk seq -L 500 -N | \
	seqkit seq head -n 25000 > data/hg38_long_mrna.fa

data/yeast.fa:
	curl http://ftp.ensemblgenomes.org/pub/fungi/release-60/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz | gzip -cdf > data/yeast.fa

data/e_coli.fa:
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gzip -cdf >data/e_coli.fa

data/soybean.fa:
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz | gzip -cdf >data/soybean.fa

data/lambda_phage.fa:
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/GCF_000840245.1_ViralProj14204_genomic.fna.gz | gzip -cdf >data/lambda_phage.fa
