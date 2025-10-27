.PHONY: all
all: dockerfiles
	# bash ./pack-alpine-tgz.sh
	bash ./pack-deb.sh

.PHONY: update-containers
update-containers:
	$(MAKE) -C dockerfiles clean
	$(MAKE) -C dockerfiles all

.PHONY: dockerfiles
dockerfiles:
	$(MAKE) -C dockerfiles all
