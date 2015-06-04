#! /usr/bin/env python
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

# I don't use the ez_setup module because it causes us to automatically build
# and install a new setuptools module, which I'm not interested in doing.

from setuptools import setup

setup (
    name = 'scanalyzer',
    version = '1.2',

    # This package actually *is* zip-safe, but I've run into issues with
    # installing it as a Zip: in particular, the install sometimes fails with
    # "bad local file header", and reloading a module after a reinstall in
    # IPython gives an ImportError with the same message. These are annoying
    # enough and I don't really care so we just install it as flat files.
    zip_safe = False,

    packages = [
        'scanalyzer',
    ],

    package_data = {
        'scanalyzer': ['scanalyzer.ui'],
    },

    install_requires = [
        'numpy >= 1.6',
    ],

    entry_points = {
        'console_scripts': [
            'scanalyzer = scanalyzer.cli:commandline',
        ],
    },

    author = 'Peter Williams',
    author_email = 'peter@newton.cx',
    description = 'Foo',
    license = 'MIT',
    keywords = 'astronomy science',
    url = 'https://newton.cx/~peter/',

    long_description = \
    '''The scanalyzer.
    ''',

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)
