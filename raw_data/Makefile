all: \
	ce11.mRNA.fa \
	ce11_chr1.fa \
	small_test.fa.fai \
	lambda_phage.fa \

small_test.fa.fai: small_test.fa
	samtools faidx small_test.fa

ce11.mRNA.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/refMrna.fa.gz | gunzip -cdf > ce11.mRNA.fa

ce11_chr1.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrI.fa.gz  | gunzip -cdf > ce11_chr1.fa
	samtools faidx ce11_chr1.fa

lambda_phage.fa:
	curl -L https://db.cngb.org/sequence/down_sequence_txt/NC_001416.1/ > lambda_phage.fa
	samtools faidx lambda_phage.fa
