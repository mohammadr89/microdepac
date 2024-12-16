# Code can be compiled using 'python3 setup.py build_ext --inplace'

build_folder = '../build'

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('radiation',
                             sources=['radiation.pyx', '../src_test/Radiation_solver.cpp'],
                             language='c++',
                             extra_compile_args=['-O3', '-std=c++14', '-march=native', '-DBOOL_TYPE=signed char', '-fno-wrapv'],
                             include_dirs=['../include', '../include_test', numpy.get_include()],
                             library_dirs=['/usr/local/Cellar/gcc/9.3.0_1/lib/gcc/9/'],
                             libraries=['gfortran', 'netcdf'],
                             extra_objects=[
                                 '{}/src/librte_rrtmgp.a'.format(build_folder),
                                 '{}/src_fortran/librte_rrtmgp_kernels.a'.format(build_folder)] )]
)
