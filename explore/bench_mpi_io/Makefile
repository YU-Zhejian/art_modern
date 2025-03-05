JOBS ?= 20

.PHONY: build
build:
	mkdir -p opt/build
	env -C opt/build cmake $(CURDIR) \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		$(CMAKE_FLAGS)
	cmake --build opt/build -j$(JOBS)
