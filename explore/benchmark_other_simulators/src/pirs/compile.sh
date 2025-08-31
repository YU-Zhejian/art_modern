g++ -w -mtune=native -O3 -msse2 -fopenmp -DSFMT_MEXP=19937 -DHAVE_CONFIG_H -DPKGDATADIR='"/usr/local/share/pirs"' -ISFMT-src-1.4 -o pirs BaseCallingProfile.cpp GCBiasProfile.cpp IndelProfile.cpp InputStream.cpp main.cpp MaskQvalsByEamss.cpp OutputStream.cpp pirs_diploid.cpp pirs_simulate.cpp Profile.cpp util.cpp SFMT-src-1.4/SFMT.c -lz -lpthread
PROFILES=""

./pirs simulate -A dist test/ref_seq.fa -M lowercase -m 500 -l 100 -x 10 -v 40 -t 20 \
    -B Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz \
    --no-indels \
    --no-gc-bias \
    -o Illumina \
    -c text
