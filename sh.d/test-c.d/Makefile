.SUFFIXES: # Remove default suffixes

CFLAGS ?= \
	-O3 -mtune=native -march=native \
	-Wall -Wextra -Wpedantic \
	-std=gnu11 \
	$(shell pkg-config --cflags zlib libxxhash libdeflate) \
	-isystem src/
CC ?= gcc

LIBS ?= $(shell pkg-config --libs zlib libxxhash libdeflate) -lm -lc

.PHONY: bin
bin: bin/fast-classify bin/find-dups bin/generate_many_contigs bin/generate_large_contigs

src/fast-classify: src/fast-classify.o
	$(CC) -o $@ $^ $(LIBS)

src/find-dups: src/find-dups.o src/bloom/bloom.o
	$(CC) -o $@ $^ $(LIBS)

src/generate_many_contigs: src/generate_many_contigs.o
	$(CC) -o $@ $^ $(LIBS)

src/generate_large_contigs: src/generate_large_contigs.o
	$(CC) -o $@ $^ $(LIBS)

large_contigs.fa.gz: src/generate_large_contigs
	src/generate_large_contigs | bgzip --threads 20 --force --index -l 9 -o large_contigs.fa.gz

large_contigs.fa.gz.fai: large_contigs.fa.gz
	samtools faidx large_contigs.fa.gz

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

bin/%: src/%
	install -D -m 755 $< $@

