
.PHONY: build
build:
	mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -G Ninja .. && ninja -j40 && cd ..

.PHONY: fmt
fmt:
	bash fmt.sh

.PHONY: touch
touch:
	bash touch-all.sh

.PHONY: testsmall
testsmall: build
	bash test_small.sh

.PHONY: testse
testse: build raw_data
	build/art \
		--qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
		--qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
		--seq_file raw_data/MANE.GRCh38.v1.3.refseq_rna.fna \
		--out_file_prefix ./tmp/single_end_com \
		--no_sam \
		--read_len 125 \
		--fcov 1

.PHONY: testpe
testpe: raw_data build
	build/art \
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
