'''
Setup script
'''

import setuptools
from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize
import numpy as np
import os

# Read requirements from requirements.txt
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Specify the package names you want to include
packages = find_packages(where="src", include=["Lyman-alpha-feedback*"])
                         
setup(
    name='anaxagoras',
    version='0.0.1',
    author='Authors',
    author_email='corresponding.author@somemail.com',
    packages=packages,
    package_dir={"": "src"},
    package_data={
        '': ['input_data/**/*'],  # Include all files in all subfolders of input_data
    },
    install_requires=requirements,
    include_package_data=True,
    include_dirs=[np.get_include()],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
