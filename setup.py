# This setup is borrowed from https://github.com/pybind/cmake_example

import os
import re
import sys
import subprocess

from setuptools import (
    setup,
    find_packages,
    Extension,
)

from setuptools.command.build_ext import build_ext


# Version number
MAJOR = 0
MINOR = 1


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
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DPYTHON_EXECUTABLE={}".format(sys.executable)
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]
        env = os.environ.copy()
        env["CXXFLAGSS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""),
            self.distribution.get_version()
        )

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)


setup(
    name="SVMTK",
    version="{0}.{1}".format(MAJOR, MINOR),
    description="A collection of tools for volume and surface meshing",
    long_description="",
    ext_modules=[CMakeExtension("brainmesh")],
    cmdclass=dict(build_ext=CMakeBuild),
    packages=["source"],
    # package_dir={"": "source"},
    entry_points={
        "console_scripts": [
            "brainmesh-convert = source.convert:main",
        ]
    },
    zip_safe=False
)
