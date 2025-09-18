DOCKER ?= docker
SINGULARITY ?= singularity

.PHONY: all clean

all: almalinux-8.sif \
	 almalinux-9.sif \
	 debian-12.sif \
	 debian-13.sif \
	 ubuntu-2404.sif \

clean:
	rm -f *.sif
	rm -f *.tar
	rm -f *.timestamp

# Generate Singularity build image from Dockerfiles

%.timestamp: %.Dockerfile
	$(DOCKER) build -t $(basename $<):latest -f $< .
	touch $@

%-squashed.tar: %.timestamp
	docker-squash $(basename $<):latest \
		--tag $(basename $<) \
		--output-path $@ \
		--load-image False
	# Remove the intermediate image
	$(DOCKER) rmi $(basename $<):latest

%.sif: %-squashed.tar
	$(SINGULARITY) build --force $@ docker-archive://$^
