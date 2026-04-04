.SUFFIXES: # Remove default suffixes

PKGCONF ?= pkgconf

# Test whether MKL provides mkl-sdl.pc or mkl-sdl-lp64.pc
MKL_SDL_PC := $(shell \
	if $(PKGCONF) --exists mkl-sdl; then \
		echo mkl-sdl; \
	elif $(PKGCONF) --exists mkl-sdl-lp64; then \
		echo mkl-sdl-lp64; \
	fi)

ifeq ($(strip $(MKL_SDL_PC)),)
$(error pkg-config could not find mkl-sdl.pc or mkl-sdl-lp64.pc)
endif

CFLAGS ?= \
	-O3 -mtune=native -march=native \
	-Wall -Wextra -Wpedantic \
	-std=gnu11 \
	$(shell $(PKGCONF) --cflags zlib libxxhash libdeflate $(MKL_SDL_PC)) \
	-isystem src/
CC ?= gcc

LIBS ?= $(shell $(PKGCONF) --libs zlib libxxhash libdeflate $(MKL_SDL_PC)) -lm -lc -fopenmp

.PHONY: all
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

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

bin/%: src/%
	install -D -m 755 $< $@

