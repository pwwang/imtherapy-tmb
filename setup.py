
# -*- coding: utf-8 -*-

# DO NOT EDIT THIS FILE!
# This file has been autogenerated by dephell <3
# https://github.com/dephell/dephell

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = ''

setup(
    long_description=readme,
    name='imtherapy-mut',
    version='0.0.0',
    description='Mutation related feature transformation module for imtherapy',
    python_requires='==3.*,>=3.7.0',
    author='pwwang',
    author_email='pwwang@pwwang.com',
    license='MIT',
    entry_points={"imtherapy_feature_transform": ["mut = imtherapy_mut:FeatureTransformMut"]},
    packages=['imtherapy_mut'],
    package_dir={"": "."},
    package_data={},
    install_requires=['imtherapy', 'pipen'],
)
