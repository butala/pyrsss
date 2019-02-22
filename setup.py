#!/usr/bin/env python

import os
from setuptools import setup, Extension, find_packages

from Cython.Build import cythonize


if 'GPSTK_SRC' in os.environ:
    assert 'GPSTK_BUILD' in os.environ
    core_lib_path = os.path.join(os.environ['GPSTK_SRC'],
                                 'core',
                                 'lib')
    gpstk_ext = Extension('pyrsss.gpstk',
                          [os.path.join('pyrsss', 'gpstk.pyx')],
                          include_dirs=[os.path.join(core_lib_path,
                                                     'GNSSCore'),
                                        os.path.join(core_lib_path,
                                                     'Utilities'),
                                        os.path.join(core_lib_path,
                                                     'Math'),
                                        os.path.join(core_lib_path,
                                                     'Math',
                                                     'Vector'),
                                        os.path.join(core_lib_path,
                                                     'RefTime'),
                                        os.path.join(core_lib_path,
                                                     'GNSSEph'),
                                        os.path.join(core_lib_path,
                                                     'TimeHandling'),
                                        os.path.join(core_lib_path,
                                                     'FileHandling'),
                                        os.path.join(core_lib_path,
                                                     'FileHandling',
                                                     'RINEX'),
                                        os.path.join(core_lib_path,
                                                     'FileHandling',
                                                     'RINEX3')],
                          library_dirs=[os.environ['GPSTK_BUILD']],
                          runtime_library_dirs=[os.environ['GPSTK_BUILD']],
                          libraries=['gpstk'],
                          language='c++')
    ext_modules = cythonize([gpstk_ext])
else:
    ext_modules = []


setup(name='pyrsss',
      version='0.2',
      description='Remote sensing and space science python tools.',
      author='Mark D. Butala, Matthew A. Grawe, et al.',
      author_email='butala@illinois.edu',
      packages=find_packages(),
      ext_modules=ext_modules)
