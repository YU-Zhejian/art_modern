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
	ce11.mRNA_head.cov_stranded.tsv \
	ce11.mRNA_head.cov_strandless.tsv \
	ce11.mRNA_head.pbsim3.transcript \


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
	seqkit grep -s -v -p N ce11.mRNA.fa | seqkit seq -m 200 | seqtk sample /dev/stdin 1000 > ce11.mRNA_head.fa

# Grouped target: one gen_cov.py run writes all three. Generated here rather than on demand from
# test-small.sh, where up to 5 concurrent test-build.py workers raced to create the same files.
ce11.mRNA_head.cov_stranded.tsv ce11.mRNA_head.cov_strandless.tsv ce11.mRNA_head.pbsim3.transcript &: \
		ce11.mRNA_head.fa
	python ../../sh.d/test-small.sh.d/gen_cov.py ce11.mRNA_head 5

ce11_chr1.fa:
	curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrI.fa.gz  | gunzip -cdf > ce11_chr1.fa

lambda_phage.fa:
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/GCF_000840245.1_ViralProj14204_genomic.fna.gz | gunzip -cdf > lambda_phage.fa
