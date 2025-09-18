FROM debian:12-slim

# Install necessary packages
RUN apt-get update && \
    apt-get full-upgrade -y && \
    apt-get install --no-install-recommends -y \
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
    devscripts \
    lintian \
    dh-cmake \
    libhts-dev \
    pkgconf \
    liblzma-dev \
    libbz2-dev \
    libfmt-dev \
    libconcurrentqueue-dev \
    libabsl-dev \
    lowdown \
    && apt clean

COPY build_deb.sh /
RUN mkdir -p /mnt/debian /mnt/build_deb
