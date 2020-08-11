import os
import re
import sys
import subprocess

from setuptools import setup,find_packages,Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand
from shutil import copyfile, copymode

class CatchTestCommand(TestCommand):
      def distutils_dir_name(self, dname):
          dir_name = "{dirname}.{platform}-{version[0]}.{version[1]}"
          return dir_name.format(dirname=dname,
                     platform=sysconfig.get_platform(),
                     version=sys.version_info)
      def run(self):
         # Run Python tests
         super(CatchTestCommand, self).run()
         print("\nPython tests complete, now running C++ tests...\n")
         # Run catch tests
         subprocess.call(['./*_test'],
                    cwd=os.path.join('build',self.distutils_dir_name('temp')),shell=True)


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
        test_bin = os.path.join(self.build_temp, 'SVMTK_test')
        self.copy_test_file(test_bin)
        print()

    def copy_test_file(self, src_file):
        '''
        Copy ``src_file`` to `tests/bin` directory, ensuring parent directory 
        exists. Messages like `creating directory /path/to/package` and
        `copying directory /src/path/to/package -> path/to/package` are 
        displayed on standard output. Adapted from scikit-build.
        '''
        # Create directory if needed
        dest_dir = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), 'tests', 'bin')
        if dest_dir != "" and not os.path.exists(dest_dir):
            print("creating directory {}".format(dest_dir))
            os.makedirs(dest_dir)

        # Copy file
        dest_file = os.path.join(dest_dir, os.path.basename(src_file))
        print("copying {} -> {}".format(src_file, dest_file))
        copyfile(src_file, dest_file)
        copymode(src_file, dest_file)

setup(
    name="SVMTK",
    version="{0}.{1}".format(MAJOR, MINOR),
    description="A collection of tools for volume and surface meshing",
    long_description="",
    ext_modules=[CMakeExtension("SVMTK")],
    cmdclass=dict(build_ext=CMakeBuild),
    test_suite='tests',
    zip_safe=False
)

