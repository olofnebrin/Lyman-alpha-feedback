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
                         
setup(
    name='Lyman_alpha_feedback',
    version='1.0.0',
    author='Authors',
    author_email='corresponding.author@somemail.com',
    packages=find_packages(where="src"),
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
