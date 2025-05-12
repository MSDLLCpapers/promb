#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os
from io import open

about = {}
# Read version number from __version__.py (see PEP 396)
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'promb', '__version__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

# Read contents of readme file into string
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def read_requirements():
    with open(os.path.join(here, 'requirements.txt'), encoding='utf-8') as f:
        return f.read().splitlines()


setup(
    name='promb',
    version=about['__version__'],
    description='promb: protein mutation burden',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='David Prihoda',
    packages=find_packages(include=['promb.*', 'promb']),
    author_email='david.prihoda@gmail.com',
    license='MIT',
    python_requires=">=3.8",
    keywords='antibody humanness, de-novo humanness, human peptide content, oasis, biophi, mutation burden, promb',
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    entry_points={
        'console_scripts': [
            'promb = promb.cli:main'
        ]
    },
    install_requires=read_requirements(),
    include_package_data=True,
    url='https://github.com/MSDLLCPapers/promb'
)
