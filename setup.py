#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup
from Cython.Build import cythonize

with open("README.rst", 'r') as f:
    long_description = f.read()

setup(
   name='opentiva',
   version='1.0',
   description='Simulation of target controlled infusions using'
               'pharmacokinetic and pharmacodynamic models.',
   author='James Selby',
   author_email='opentiva@jpselby.co.uk',
   project_urls={
        'Documentation': 'https://opentiva.readthedocs.io',
        'Source': 'https://github.com/opentiva/opentiva',
   },
   license='LICENSE',
   packages=['opentiva'],
   ext_modules=cythonize('opentiva/pkpd.pyx'),
   zip_safe=False,
   install_requires=['numpy', 'scipy', 'cython'],
   classifiers=[
       "Development Status :: 5 - Production/Stable",
       "Intended Audience :: Healthcare Industry",
       "Intended Audience :: Science/Research",
       "License :: OSI Approved :: MIT License",
       "Natural Language :: English",
       "Operating System :: OS Independent",
       "Programming Language :: Python :: 3",
       "Programming Language :: Python :: 3.8",
       "Programming Language :: Python :: 3.9",
       "Programming Language :: Python :: Implementation :: CPython",
       "Programming Language :: Python :: Implementation :: PyPy"
       ],
)
