FROM debian:buster

WORKDIR /svmtk

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
    python3-pip

# Downlod CGAL and pybind
RUN mkdir external && cd external && \
    git clone https://github.com/pybind/pybind11.git --branch=v2.6.2 && \
    git clone https://github.com/CGAL/cgal.git --branch=v5.3.0

RUN python3 -m pip install setuptools

RUN python3 -m pip install numpy

# ADD demos demos
ADD examples examples
ADD include include
ADD python python
# ADD source source
ADD tests tests
ADD CMakeLists.txt .
ADD setup.py .
ADD README.md .
ADD REQUIREMENTS.md .

RUN python3 setup.py install
ADD tests tests 

