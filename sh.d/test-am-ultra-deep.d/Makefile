CFLAGS ?= -O3 -mtune=native -march=native -Wall -Wextra -Wpedantic
CC ?= gcc

LIBS ?= -lz -ldeflate

.PHONY: all
all: fast-classify

fast-classify: fast-classify.o getdelim.o
	$(CC) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm *.o fast-classify
