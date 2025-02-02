CMAKE_FLAGS ?= 

.PHONY: build
build:
	mkdir -p opt/build_debug
	env -C opt/build_debug cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build opt/build_debug -j40
	env -C opt/build_debug ctest --output-on-failure
	opt/build_debug/art_modern --help
	opt/build_debug/art_modern --version # mpiexec --verbose -n 5 

.PHONY: release
release:
	mkdir -p opt/build_release
	env -C opt/build_release cmake \
		-DCMAKE_BUILD_TYPE=Release \
		-DCEU_CM_SHOULD_USE_NATIVE=ON \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build opt/build_release -j40
	# cpack --config opt/build_release/CPackSourceConfig.cmake

.PHONY: rel_with_dbg_alpine
rel_with_dbg_alpine:
	mkdir -p opt/build_rel_with_dbg_alpine
	env -C opt/build_rel_with_dbg_alpine cmake \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCEU_CM_SHOULD_ENABLE_TEST=OFF \
		-DCEU_CM_SHOULD_USE_NATIVE=OFF \
		-DBUILD_SHARED_LIBS=OFF \
		-DUSE_RANDOM_GENERATOR=BOOST \
        $(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build opt/build_rel_with_dbg_alpine -j40

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

.PHONY: testsmall-release
testsmall-release: release raw_data
	env ART=opt/build_release/art_modern bash sh.d/test_small.sh

.PHONY: raw_data
raw_data:
	$(MAKE) -C data/raw_data

.PHONY: clean
clean:
	rm -fr opt tmp build
