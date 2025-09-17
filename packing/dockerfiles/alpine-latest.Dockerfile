FROM alpine:latest

# Install necessary packages
RUN apk update && apk upgrade --no-cache && \
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
    htslib-static \
    htslib-dev \
    xz-static \
    xz-dev \
    bzip2-static \
    bzip2-dev \
    libdeflate-static \
    libdeflate-dev
# Bundled {fmt}, moodycamel, and abseil libraries are used \
