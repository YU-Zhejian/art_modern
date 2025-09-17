FROM debian:12-slim

COPY build_deb.sh /build_deb.sh

# Install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    binutils \
    libboost-all-dev \
    zlib1g-dev \
    make \
    python3 \
    sed \
    grep \
    cmake \
    dh-cmake \
    libhts-dev \
    pkgconf \
    liblzma-dev \
    libbz2-dev \
    libfmt-dev \
    libconcurrentqueue-dev \
    libabsl-dev \
    && apt clean \
