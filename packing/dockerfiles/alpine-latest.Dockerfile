FROM alpine:latest

# Install necessary packages
RUN apk update && \
    apk upgrade --no-cache && \
    apk add --no-cache \
    g++ \
    binutils \
    boost-dev \
    zlib-dev \
    make \
    python3 \
    cmake \
    zlib-static \
    icu-static \
    coreutils \
    boost1.84-static \
    sed \
    grep \
    xz-static \
    xz-dev \
    bzip2-static \
    bzip2-dev \
    libdeflate-static \
    libdeflate-dev
# Bundled {fmt}, moodycamel, and HTSLib libraries are used
# Since we're not releasing a package, using bundled libraries is acceptable.

COPY build_alpine_tgz.sh /
RUN mkdir -p /mnt/art_modern-master
