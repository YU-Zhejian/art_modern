all: \
	MANE.GRCh38.v1.3.refseq_rna.fna \
	ce11_chr1.fa \
	MANE.GRCh38.v1.3.refseq_genomic.gtf \
	MANE.GRCh38.v1.3.summary.txt \


MANE.GRCh38.v1.3.refseq_rna.fna:
	curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.refseq_rna.fna.gz | gunzip -cdf > MANE.GRCh38.v1.3.refseq_rna.fna
	# salmon index -t MANE.GRCh38.v1.3.refseq_rna.fna -i MANE.GRCh38.v1.3.refseq_rna.salmon_idx -p 40

ce11_chr1.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrI.fa.gz  | gunzip -cdf > ce11_chr1.fa

MANE.GRCh38.v1.3.refseq_genomic.gtf:
	curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.refseq_genomic.gtf.gz | gunzip -cdf > MANE.GRCh38.v1.3.refseq_genomic.gtf

MANE.GRCh38.v1.3.summary.txt:
	curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.summary.txt.gz | gunzip > MANE.GRCh38.v1.3.summary.txt
	# cat MANE.GRCh38.v1.3.summary.txt | cut -f 4,6 | sed -E 's;(\S+)\s(\S+);\2\t\1;' | sed '1d' > MANE_salmon_genemap.tsv