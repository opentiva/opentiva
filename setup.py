#!/usr/bin/env python
# encoding: utf-8

from setuptools import Extension, setup

try:
    from Cython.Build import cythonize
except ImportError:
    use_cython = False
else:
    use_cython = True

ext = '.pyx' if use_cython else '.c'

extensions = [Extension("opentiva.pkpd", ["opentiva/pkpd"+ext])]

if use_cython:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

with open("README.rst", 'r') as f:
    readme = f.read()

setup(
   name='opentiva',
   version='1.0.4',
   description='Simulation of target controlled infusions using'
               'pharmacokinetic and pharmacodynamic models.',
   long_description=readme,
   long_description_content_type='text/x-rst',
   author='James Selby',
   author_email='opentiva@jpselby.co.uk',
   project_urls={
        'Documentation': 'https://opentiva.readthedocs.io',
        'Source': 'https://github.com/opentiva/opentiva',
   },
   license='LICENSE',
   packages=['opentiva'],
   ext_modules=extensions,
   zip_safe=False,
   include_package_data=True,
   install_requires=['numpy', 'scipy'],
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
