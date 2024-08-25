all: \
	ce11.mRNA.fa \
	ce11.mRNA.fa.fai \
	ce11_chr1.fa \
	ce11_chr1.fa.fai \
	small_test.fa.fai \
	lambda_phage.fa \
	lambda_phage.fa.fai \
	ce11.mRNA_head.fa \
	ce11.mRNA_head.fa.fai \


ce11.mRNA.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/refMrna.fa.gz | gunzip -cdf > ce11.mRNA.fa
ce11.mRNA_head.fa:
	seqtk sample ce11.mRNA.fa 1000 > ce11.mRNA_head.fa

ce11_chr1.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrI.fa.gz  | gunzip -cdf > ce11_chr1.fa

lambda_phage.fa:
	curl -L https://db.cngb.org/sequence/down_sequence_txt/NC_001416.1/ > lambda_phage.fa

%.fa.fai: %.fa
	samtools faidx $<