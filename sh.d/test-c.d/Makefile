CFLAGS ?= -O3 -mtune=native -march=native -Wall -Wextra -Wpedantic  -std=c11
CC ?= gcc

.PHONY: all
all: src/generate_many_contigs

large_contigs.fa.gz: src/generate_large_contigs
	src/generate_large_contigs | bgzip --threads 20 --force --index -l 9 -o large_contigs.fa.gz

large_contigs.fa.gz.fai: large_contigs.fa.gz
	samtools faidx large_contigs.fa.gz

%: %.c
	$(CC) $(CFLAGS) -o $@ $<
