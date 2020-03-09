FROM debian:buster

WORKDIR /svmtk

RUN apt update && apt install -y \
    bzip2 \
    cmake \
    curl \
    g++ \
    libboost-all-dev \
    libeigen3-dev \
    libgmp3-dev \
    libmpfr-dev \
    xz-utils \
    zlib1g-dev \
    git

RUN curl -sL https://github.com/CGAL/cgal/releases/download/releases/CGAL-4.13/CGAL-4.13.tar.xz | tar xpvfJ -
RUN cd CGAL-4.13 && cmake -DWITH_Eigen3:BOOL=ON . && make && cat CMakeCache.txt
RUN git clone --recursive https://github.com/SVMTK/SVMTK.git --branch=docker
