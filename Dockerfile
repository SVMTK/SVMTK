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
RUN curl -sL https://github.com/CGAL/cgal/releases/download/releases/CGAL-4.13/CGAL-4.13.tar.xz | tar xpvfJ -
RUN cd CGAL-4.13 && cmake -DWITH_Eigen3:BOOL=ON . && make && make install

RUN mkdir external && cd external && \
    git clone https://github.com/pybind/pybind11.git --branch=v2.4.3 && \
    cd $HOME

RUN python3 -m pip install setuptools

ADD demos demos
ADD examples examples
ADD include include
ADD python python
ADD source source
ADD test test
ADD CMakeLists.txt .
ADD setup.py .

RUN python3 setup.py install
