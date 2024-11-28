.PHONY: build
build:
	mkdir -p build
	env -C build cmake \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		-G Ninja ..
	env -C build ninja -j40
	env -C build ctest --output-on-failure
	build/art_modern --help
	mpiexec --verbose -n 5 build/art_modern --version

.PHONY: fmt
fmt:
	bash sh.d/fmt.sh

.PHONY: scc
scc:
	bash sh.d/scc.sh

.PHONY: touch
touch:
	bash sh.d/touch-all.sh

.PHONY: testsmall
testsmall: build raw_data
	bash sh.d/test_small.sh

.PHONY: raw_data
raw_data:
	$(MAKE) -C raw_data

.PHONY: clean
clean:
	rm -fr build build_profile tmp
