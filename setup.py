#!/usr/bin/env python

from distutils.core import setup

setup(name='pyrsss',
      version='0.1',
      description='Remote sensing and space science python tools.',
      author='Mark D. Butala, Matthew A. Grawe, et al.',
      author_email='butala@illinois.edu',
      package_dir = {'pyrsss': '.'},
      packages=['pyrsss',
                'pyrsss.emission',
                'pyrsss.gps',
                'pyrss.ionex',
                'pyrsss.mag',
                'pyrsss.util'])
