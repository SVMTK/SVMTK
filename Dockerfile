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
    git \
    python3-pip

# Istall and compile CGAL
ARG CGAL_VERSION
ENV CGAL_VERSION 5.0.2
RUN curl -sL https://github.com/CGAL/cgal/releases/download/releases/CGAL-${CGAL_VERSION}/CGAL-${CGAL_VERSION}.tar.xz | tar xpvfJ -
RUN cd CGAL-${CGAL_VERSION} && cmake -DWITH_Eigen3:BOOL=ON . && make && make install

RUN mkdir external && cd external && \
    git clone https://github.com/pybind/pybind11.git --branch=v2.4.3 && \
    git clone https://github.com/catchorg/Catch2.git

RUN python3 -m pip install setuptools

# ADD demos demos
ADD examples examples
ADD include include
ADD python python
# ADD source source
ADD test test
ADD CMakeLists.txt .
ADD setup.py .

RUN python3 setup.py install
