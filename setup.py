#!/usr/bin/env python

from setuptools import setup

setup(
  name='BTM',
  version='1.0.15',

  author='Shuzhao Li',
  author_email='shuzhao.li@gmail.com',
  description='Blood Transcription Modules for transcriptomics analysis',
  long_description=open('README.md').read(),
  url='https://github.com/shuzhao-li/BTM',
  license='MIT',

  keywords='transcriptomics analysis bioinformatics immunology systems biology',

  classifiers=[
    'Development Status :: 5 - Production/Stable',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=['BTM'],
  include_package_data=True,
  zip_safe=True,

  install_requires=[
    'numpy',
    'scipy',
  ],

)
