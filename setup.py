import platform
import sys

from setuptools import (
    setup,
    find_packages,
)

# Version number
major = 0
minor = 1


setup(
    name="brainmesh",
    version="{0}.{1}".format(major, minor),
    description="A collection of scripts for brain meshing.",
    packages=find_packages("source"),
    package_dir={"": "source"},
    entry_points={
        "console_scripts": [
            "asc2domain = brainmesh.utils.asc2domain:main",
        ]
    }
)
