import os
from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(file_name):
    return open(os.path.join(os.path.dirname(__file__), file_name)).read()


with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="pvt_tools",
    version="0.0.4",
    author="Equinor ASA",
    author_email="pih@equinor.com",
    description="Useful PVT functions in python",
    long_description=open("README.md").read(),
    packages=["pvt_tools"],
    package_dir={"": "src"},
    test_suite="tests",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
    ],
    install_requires=requirements,
)
