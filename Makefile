.PHONY: build
build:
	mkdir -p build
	env -C build cmake -DCMAKE_BUILD_TYPE=Release -G Ninja ..
	env -C build ninja -j40
	env -C build ctest

.PHONY: build_external_htslib
build_external_htslib:
	mkdir -p build_external_htslib
	env -C build_external_htslib cmake \
		-DCMAKE_BUILD_TYPE=Release \
		-DUSE_HTSLIB="hts" \
		-G Ninja ..
	env -C build_external_htslib ninja -j40
	env -C build_external_htslib ctest

.PHONY: fmt
fmt:
	bash fmt.sh

.PHONY: scc
scc:
	bash scc.sh

.PHONY: touch
touch:
	bash touch-all.sh

.PHONY: testsmall
testsmall: build raw_data
	bash test_small.sh

.PHONY: testse
testse: build raw_data
	build/art_modern \
		--qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
		--qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
		--seq_file raw_data/MANE.GRCh38.v1.3.refseq_rna.fna \
		--out_file_prefix ./tmp/single_end_com \
		--no_sam \
		--read_len 125 \
		--fcov 1

.PHONY: testpe
testpe: raw_data build
	build/art_modern \
		--qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
		--qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
		--seq_file raw_data/ce11_chr1.fa \
		--out_file_prefix ./tmp/pe_genomic \
		--read_len 125 \
		--fcov 5.0 \
		--pe_frag_dist_mean 200 \
		--pe_frag_dist_std_dev 10 \
		--is_pe \
		--parallel_on_read


.PHONY: raw_data
raw_data:
	$(MAKE) -C raw_data

.PHONY: build-alpine
build-alpine:
	cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF; ninja

.PHONY: sync-ceu-cm
sync-ceu-cm:
	rm -fr deps/cmake_collections
	cp -r ../libceu/cmake_collections/ deps/
