.PHONY: all
all: dockerfiles
	bash ./pack-deb.sh
	bash ./pack-alpine_tgz.sh

.PHONY: update-containers
update-containers:
	$(MAKE) -C dockerfiles clean
	$(MAKE) -C dockerfiles all

.PHONY: dockerfiles
dockerfiles:
	$(MAKE) -C dockerfiles all
