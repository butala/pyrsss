#!/usr/bin/env python

from glob import glob

from numpy.distutils.core import setup, Extension

madrigal_ext = Extension(name='madrigal',
                         sources=glob('madrigal/*.f') + \
                                 ['madrigal/madrigal.pyf'])


setup(name='pyrsss',
      version='0.1',
      description='Remote sensing and space science python tools.',
      author='Mark D. Butala, Matthew A. Grawe, et al.',
      author_email='butala@illinois.edu',
      package_dir = {'pyrsss': '.'},
      packages=['pyrsss',
                'pyrsss.emission',
                'pyrsss.gps',
                'pyrsss.ionex',
                'pyrsss.mag',
                'pyrsss.util'],
      ext_modules=[madrigal_ext])
