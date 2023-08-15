#!/usr/bin/env python

import os
from setuptools import setup, Extension, find_packages

from Cython.Build import cythonize


if 'GNSSTK_SRC' in os.environ:
    assert 'GNSSTK_BUILD' in os.environ
    core_lib_path = os.path.join(os.environ['GNSSTK_SRC'],
                                 'core',
                                 'lib')
    gnsstk_ext = Extension('pyrsss.gnsstk',
                          [os.path.join('pyrsss', 'gnsstk.pyx')],
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
                                                      'RINEX3'),
                                         os.environ['GNSSTK_BUILD']],
                           extra_compile_args=['-std=c++11'],
                           library_dirs=[os.environ['GNSSTK_BUILD']],
                           runtime_library_dirs=[os.environ['GNSSTK_BUILD']],
                           libraries=['gnsstk'],
                           language='c++')
    ext_modules = cythonize([gnsstk_ext])
else:
    ext_modules = []


setup(name='pyrsss',
      version='0.2',
      description='Remote sensing and space science python tools.',
      author='Mark D. Butala, Matthew A. Grawe, et al.',
      author_email='butala@illinois.edu',
      packages=find_packages(),
      ext_modules=ext_modules)
