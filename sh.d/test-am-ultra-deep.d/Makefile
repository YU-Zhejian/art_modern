CFLAGS ?= -O3 -mtune=native -march=native -Wall -Wextra -Wpedantic -g
CC ?= gcc

LIBS ?= -lz -ldeflate -lxxhash

.PHONY: all
all: fast-classify find-dups

fast-classify: fast-classify.o netbsd/getdelim.o
	$(CC) -o $@ $^ $(LIBS)

find-dups: find-dups.o netbsd/getdelim.o bloom/bloom.o
	$(CC) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f *.o */*.o fast-classify find-dups
