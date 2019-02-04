FROM ubuntu:18.04

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
    zlib1g-dev \
    git \
    python3-setuptools

RUN curl -sL https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13/CGAL-4.13.tar.xz | tar -xJf -
RUN cd CGAL-4.13 && cmake -DWITH_Eigen3:BOOL=ON . && make && make install
RUN git clone --recursive https://github.com/jakobes/BrainMesh.git && cd BrainMesh && python3 setup.py install
