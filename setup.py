import os
import re
import sys
import subprocess

from setuptools import setup,find_packages,Extension
from setuptools.command.build_ext import build_ext 
from setuptools.command.test import test
from shutil import copyfile, copymode

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


# Version number
MAJOR = 1
MINOR = 0

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            names = ", ".join(e.name for e in self.extensions)
            msg = "CMake must be installed to build the following extensions: {}".format(names)
            raise RuntimeError(msg)

        for ext in self.extensions:
            self.build_extension(ext)
        
        
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir)#, "-DPYTHON_EXECUTABLE={}".format(sys.executable)
        ]
        
        
        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]
        env = os.environ.copy()
   
        env["CXXFLAGSS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", "-g"),
            self.distribution.get_version()
        )
        # default to 3 build threads
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in env:
            env["CMAKE_BUILD_PARALLEL_LEVEL"] = "1"

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build","."] + build_args, cwd=self.build_temp)

     
setup(
    name="SVMTK",
    version="{}.{}".format(MAJOR, MINOR),
    description="A collection of tools for volume and surface meshing",
    url="",
    packages=find_packages(include=['SVMTK', 'SVMTK.*']),
    python_requires='>=3',
    license="GNU GENERAL PUBLIC LICENSE",
    keywords="Surface Volume Meshing Toolkit " ,
    long_description=long_description,
    long_description_content_type='text/markdown',
    ext_modules=[CMakeExtension("SVMTK")],
    cmdclass={'build_ext': CMakeBuild},
    test_suite='tests',
    zip_safe=False
)

