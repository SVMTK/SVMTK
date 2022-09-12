:orphan:

.. _Requirementsdoc:

Requirements 
============================================

The installation script of SVMTK will automatically try to find 
the necessary packages. Some of the packages can already be installed, 
so it can be practical to run the install script first to get an overview 
of the missing requirements. 

Administrator privileges can be avoided by installing in a local directory. 

Be aware that the installation may fail if multiple versions of packages exists. 

The user can also opt for the usage of docker images that was described in
the 
:ref:`Dockerdoc`.

Configure Pybind11
---------------------------

The relevant files can be found in the folder external/pybind11 after::

   git clone --recursive https://github.com/SVMTK/SVMTK

The submodule can be updated with::

   git submodule update --init

and the submodule version can be switched with::

   cd external/pybind11

   git checkout <version>

Configure CGAL
-------------------------------
The relevant files can be found in external/cgal after::

   git clone --recursive https://github.com/SVMTK/SVMTK

The submodule can be updated with:: 

   git submodule update --init

and the submodule version can be switched with:: 

   cd external/cgal

   git checkout releases/CGAL-5.3.0

Tested CGAL versions 
  - 5.3.0

Python>=3.6 
-------------------------------------------
Installing python3.6::
 
   sudo apt-get update

   sudo apt-get install -y python3.6

Check version of python3:: 

   python3 --version

In case that the version does not correspond to python3.6, we 
can change the default version python3 with::

   sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6

   sudo update-alternatives --config python3

Additional information regarding python can be read at `Python <https://www.python.org>`_.


C++>7 requriement 
--------------------------------
The compiler requirement c++1z require that gcc>=7:: 

   sudo apt-get install -y software-properties-common

   sudo add-apt-repository ppa:ubuntu-toolchain-r/test

   sudo apt-get update 

   sudo apt-get install -y g++-7

   sudo apt-get install -y build-essential

Check version of gcc in terminal::

   gcc --version

In case that the version does not correspond to gcc-7, we 
can change the default of gcc with:: 

   sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7

   sudo update-alternatives --config gcc

Boost 
---------------------------------

Installing boost using terminal::

   sudo apt-get install -y libboost-all-dev

From source 

Download and unpack a version of boost from `Download boost <https://www.boost.org/users/download/>`_.

Open the boost folder and run:: 

   ./bootstrap

Change the installation folder with:: 

   ./bootstrap.sh --prefix=new/install/folder

Complete the installation with::
 
   ./b2 install

Further details can be found at `Boost <http://www.boost.org>`_.

GMP
--------------------------------------
Installing the GMP using terminal::

   sudo apt-get install -y libgmp3-dev

From source :

Download and unpack a version of `Download GMP <https://ftp.gnu.org/gnu/gmp>`_.

Open the GMP folder and run:: 

   ./configure

   make

Change the installation folder with:: 

   ./configure --prefix=new/install/folder
   
Complete the installation with::

   make install

Further details can be found at `GMP <https://gmplib.org>`_.

MPFR
---------------------------------------
MPFR requires GMP version 5.0.0 or later.
Installing MPFR using terminal:: 

   sudo apt-get install -y libmpfr-dev 

From source :

Download and unpack a version of MPFR from `Download MPFR <https://www.mpfr.org/mpfr-current/#download>`_.

Open the MPFR folder and run:: 

   ./configure

   make

Change the installation folder with 

   ./configure --prefix=new/install/folder

Complete the installation with

   make install

Further details can be found at `MPFR <http://www.mpfr.org>`_.

CMake >=3.12
---------------------------------
Installing CMake using terminal::

   sudo apt-get install -y cmake

From source :

Download and unpack a version of CMake from `Download Cmake <http://cmake.org/download>`_.

Open the CMake folder and complete the installation with::

   ./bootstrap && make && sudo make install

Change the install folder with::

   ./bootstrap --prefix=new/install/folder && make && make install

There are more information regarding `CMake <http://cmake.org>`_.

Eigen>=3.2
------------------------------------
Installing Eigen using terminal::
 
   sudo apt-get install -y libeigen3-dev

From source :

Download and unpack a version>=3.2 from `Download Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)>`_.

Open the Eigen folder and make a build directory:: 

   mkdir build_dir

Enter the build directory::

   cd build_dir

Then use CMake to build.
 
   cmake ../

The installation folder can be changed with::
 
   cmake ../ -DCMAKE_INSTALL_PREFIX=new/install/folder

Complete the installation with:: 

   make install

Further details can be found at `Eigen <http://eigen.tuxfamily.org>`_.


.. raw:: latex

    \newpage

