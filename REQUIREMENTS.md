## Installing additional requirements 
The installation script of SVMTK will automatically try to find 
the necessary packages. Some of the packages can already be installed, 
so it can be practical to run the install script first to get an overview 
of the missing requirements. 

Administrator privileges can be avoided by installing in a local directory. 

The user can also opt for the usage of docker images that was described in
the README.md

# Python>=3.6 
Installing python3.6
 
`sudo apt-get update`

`sudo apt-get install -y python3.6`

Check version of python3 

`python3 --version`

In case that the version does not correspond to python3.6, we 
can change the default version python3 with

`sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6` 

`sudo update-alternatives --config python3`

Additional information regarding python can be read at
<https://www.python.org>

# C++>7 requriement 
The compiler requirement c++1z require that gcc>=7. 

`sudo apt-get install -y software-properties-common`

`sudo add-apt-repository ppa:ubuntu-toolchain-r/test`

`sudo apt update`

`sudo apt install -y g++-7`

`sudo apt-get install -y build-essential`

Check version of gcc in terminal

`gcc --version`

In case that the version does not correspond to gcc-7, we 
can change the default of gcc with 

`sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7`

`sudo update-alternatives --config gcc`

# Boost 
Installing boost using terminal

`sudo apt-get install -y libboost-all-dev`

From source 

Download and unpack a version of boost from
[Download](https://www.boost.org/users/download/) 

Open the boost folder and run 

`./bootstrap`

Change the installation folder with 

`./bootstrap.sh --prefix=new/install/folder`

Complete the installation with 

`./b2 install`

Further details can be found at <http://www.boost.org>

# GMP 
Installing the GMP using terminal

`sudo apt-get install -y libgmp3-dev`

From source :

Download and unpack a version of GMP from 
[Download](https://ftp.gnu.org/gnu/gmp)

Open the GMP folder and run 

`./configure`

`make`

Change the installation folder with 

`./configure --prefix=new/install/folder`

Complete the installation with

`make install`

Further details can be found at <https://gmplib.org>

# MPFR
MPFR requires GMP version 5.0.0 or later.
Installing MPFR using terminal 

`sudo apt-get install -y libmpfr-dev`

From source :

Download and unpack a version of MPFR from 
[Download](https://www.mpfr.org/mpfr-current/#download)

Open the MPFR folder and run 

`./configure`

`make`

Change the installation folder with 

`./configure --prefix=new/install/folder`

Complete the installation with

`make install`

Further details can be found at <http://www.mpfr.org>

# CMake >=3.5
Installing CMake using terminal 

`sudo apt-get install -y cmake`

From source :

Download and unpack a version of CMake from 
[Download](http://cmake.org/download)

Open the CMake folder and complete the installation with

`./bootstrap && make && sudo make install`

Change the install folder with 

`./bootstrap --prefix=new/install/folder && make && make install`

There are more information regarding CMake at <http://cmake.org>

# Eigen>=3.2
Installing Eigen using terminal

`sudo apt-get install -y libeigen3-dev`

From source :

Download and unpack a version>=3.2 from [Download](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)

Open the Eigen folder and make a build directory 

`mkdir build_dir`

Enter the build directory

`cd build_dir`

Then use CMake to build.

`cmake ../`

The installation folder can be changed with

`cmake ../ -DCMAKE_INSTALL_PREFIX=new/install/folder`

Complete the installation with 

`make install`

Further details can be found at <http://eigen.tuxfamily.org>

