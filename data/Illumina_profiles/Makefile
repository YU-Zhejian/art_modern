
# Generate a png for all txts
PNGS = $(patsubst %.txt,%.png,$(wildcard *.txt))
all: $(PNGS)


%.png: %.txt
	art-profile-fastqc -i $< -o $@ 2> $<.log --figwidth 10
