FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    bzip2 \
    cmake \
    curl \
    g++ \
    libboost-all-dev \
    libeigen3-dev \
    libgmp3-dev \
    libmpfr-dev \
    xz-utils \
    zlib1g-dev

WORKDIR foo

RUN curl -sL https://github.com/CGAL/cgal/releases/download/releases/CGAL-4.11/CGAL-4.11.tar.xz | tar xpvfJ -
RUN cd CGAL-4.11 && cmake -DWITH_Eigen3:BOOL=ON . && make && cat CMakeCache.txt
