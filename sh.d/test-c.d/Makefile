.SUFFIXES: # Remove default suffixes

CFLAGS ?= \
	-O3 -mtune=native -march=native \
	-Wall -Wextra -Wpedantic \
	-std=gnu11 \
	$(shell pkg-config --cflags zlib libxxhash libdeflate mkl-sdl-lp64) \
	-isystem src/
CC ?= gcc

LIBS ?= $(shell pkg-config --libs zlib libxxhash libdeflate mkl-sdl-lp64) -lm -lc -fopenmp

.PNONY: all
all: bin

.PHONY: bin
bin: \
	bin/fast-classify \
	bin/find-dups \
	bin/generate_many_contigs \
	bin/generate_large_contigs \
	bin/generate_rand_large_contigs

.PHONY: clean
clean:
	rm -f src/*.o src/*/*.o bin/*

src/fast-classify: src/fast-classify.o
	$(CC) -o $@ $^ $(LIBS)

src/find-dups: src/find-dups.o src/bloom/bloom.o
	$(CC) -o $@ $^ $(LIBS)

src/generate_many_contigs: src/generate_many_contigs.o
	$(CC) -o $@ $^ $(LIBS)

src/generate_large_contigs: src/generate_large_contigs.o
	$(CC) -o $@ $^ $(LIBS)

src/generate_rand_large_contigs: src/generate_rand_large_contigs.o src/twobit.o
	$(CC) -o $@ $^ $(LIBS)

data/large_contigs.fa.gz: bin/generate_large_contigs
	bin/generate_large_contigs | bgzip --threads 20 --force --index -l 9 -o data/large_contigs.fa.gz

data/large_contigs.fa.gz.fai: data/large_contigs.fa.gz
	samtools faidx data/large_contigs.fa.gz

data/rand_large_contigs.fa: bin/generate_rand_large_contigs
	bin/generate_rand_large_contigs >data/rand_large_contigs.fa

data/rand_large_contigs.fa.fai: data/rand_large_contigs.fa
	samtools faidx rand_large_contigs.fa

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

bin/%: src/%
	install -D -m 755 $< $@

