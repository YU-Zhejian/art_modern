CFLAGS ?= -O3 -mtune=native -march=native -Wall -Wextra -Wpedantic  -std=c11
CC ?= gcc

.PHONY: all
all: generate_many_contigs

%: %.c
	$(CC) $(CFLAGS) -o $@ $<
