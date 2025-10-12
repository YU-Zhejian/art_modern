.PHONY: all
all: \
	ce11.mRNA.fa \
	ce11.mRNA.fa.fai \
	ce11_chr1.fa \
	ce11_chr1.fa.fai \
	small_test.fa.fai \
	ce11.mRNA_head.fa \
	ce11.mRNA_head.fa.fai \
	ce11.mRNA_head.fa.gz \
	ce11.mRNA_head.fa.gz.fai \
	ce11_chr1.fa.gz \
	ce11_chr1.fa.gz.fai \
	lambda_phage.fa \
	lambda_phage.fa.fai \
	err.fa.fai \
	err.fa.gz \
	err.fa.gz.fai \


.SUFFIXES:

%.fa.gz: %.fa
	bgzip -i -k -f $< -c > $@

%.fa.fai: %.fa
	samtools faidx $<

%.fa.gz.fai: %.fa.gz
	samtools faidx $<

ce11.mRNA.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/refMrna.fa.gz | gunzip -cdf > ce11.mRNA.fa

ce11.mRNA_head.fa:
	seqtk sample ce11.mRNA.fa 1000 > ce11.mRNA_head.fa

ce11_chr1.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrI.fa.gz  | gunzip -cdf > ce11_chr1.fa

lambda_phage.fa:
	efetch -db nuccore -id 9626243 -format fasta > output.fa > lambda_phage.fa
