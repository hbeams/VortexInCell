#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension 
from pathlib import Path

HOME = './../../../../..' 
module1 = Extension('Hockney',
        define_macros = [('DIM', '2'),('ANIMATION', '1')],
        include_dirs = ['/usr/local',
                        './../../../../../src/wrappedCode/Hockney',
                        './../../../../../src/wrappedCode/RectMDArray',
                        './../../../../../src/wrappedCode/fftTools',
                        '/usr/include', '/usr/local/include'],
        libraries = ['python2.7', 'boost_python'],
        library_dirs = ['/usr/local/Cellar', '/usr/local/lib', '/usr/lib'],       
        sources = ['./../../../../../src/wrappedCode/fftTools/FFT1DW.cpp', './../../../../../src/wrappedCode/fftTools/FFTMD.cpp', './../../../../../src/wrappedCode/fftTools/PowerItoI.cpp','./../../../../../src/wrappedCode/Hockney/Hockney.cpp','./../../../../../src/wrappedCode/RectMDArray/DBox.cpp' ],
        extra_compile_args=['-std=c++14'])
setup(name='PackageName',
        ext_modules=[module1])
