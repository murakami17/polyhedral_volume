# coding: utf-8
# Copyright (c) 2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the MIT License.

from setuptools import setup
from setuptools import find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='polyhedral_volume',
    version='1.0.0',
    description='Calculate volume of polyhedral clusters from atomic coordinates.',
    long_description=readme,
    author='Taku MURAKAMI',
    author_email='murakami.taku.17@shizuoka.ac.jp',
    url='https://github.com/murakami17/polyhedral_volume',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

