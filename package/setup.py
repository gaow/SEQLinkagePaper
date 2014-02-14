# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import sys, os, subprocess
import imp
for item in ['faulthandler', 'numpy', 'matplotlib', 'prettyplotlib']:
    try:
        imp.find_module(item)
    except ImportError:
        sys.exit('Cannot build package: missing module "{}"!'.format(item))

from distutils.core import setup, Extension
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

from glob import glob
from source import VERSION

NAME = "SEQLinco"

if sys.platform == "win32":
    sys.exit('Windows OS is not supported.')
    
SWIG_OPTS = ['-c++', '-python', '-O', '-shadow', '-keyword',
             '-w-511', '-w-509', '-outdir', '.']

if sys.version_info.major == 2:
    PYVERSION = 'py2'
else:
    SWIG_OPTS.append('-py3')
    PYVERSION = 'py3'
#
def fixPath(fn, prefix = "source/umich"):
    if type(fn) is list:
        return [os.path.join(prefix, x) for x in fn]
    else:
        return os.path.join(prefix, fn)
#
WRAPPER_CPP = fixPath('libcore_{0}.cpp'.format(PYVERSION), "source")
WRAPPER_PY = fixPath('libcore_{0}.py'.format(PYVERSION), "source")
WRAPPER_I = fixPath('libcore.i', "source")
HEADER = fixPath(['libcore.i', 'chp.hpp', 'Core.hpp', 'VCFstream.hpp', 'Exception.hpp'], "source")
CPP = fixPath(['Core.cpp'], "source")
# generate wrapper files
try:
    ret = subprocess.call(['swig', '-python', '-external-runtime', 'swigpyrun.h'], shell=False)
    if ret != 0:
        sys.exit('Failed to generate swig runtime header file. Is "swig" installed?')
    #
    if (not os.path.isfile(WRAPPER_PY) or not os.path.isfile(WRAPPER_CPP) or \
        os.path.getmtime(WRAPPER_CPP) < max([os.path.getmtime(x) for x in HEADER + CPP])):
        ret = subprocess.call(['swig'] + SWIG_OPTS + ['-o', WRAPPER_CPP, WRAPPER_I], shell=False)
        if ret != 0:
            sys.exit('Failed to generate C++ extension.')
        os.rename('libcore.py', WRAPPER_PY)
    os.remove('swigpyrun.h')
except OSError as e:
    sys.exit('Failed to generate wrapper file: {0}'.format(e))
#
UMICH_FILES = \
  glob(fixPath("clusters/*.cpp")) + glob(fixPath("libsrc/*.cpp")) + \
  glob(fixPath("merlin/*.cpp")) +  glob(fixPath("pdf/*.cpp")) + \
  glob(fixPath("klib/*.c")) +  glob(fixPath("general/*.cpp")) + glob(fixPath("vcf/*.cpp"))
# Under linux/gcc, lib stdc++ is needed for C++ based extension.
libs = ['stdc++'] if sys.platform == 'linux2' else []
gccargs = ["-O3", "-march=native", "-std=c++11", "-D_FILE_OFFSET_BITS=64", "-D__ZLIB_AVAILABLE__ "]
#
LIBCORE_MODULE = [
    Extension('{}._libcore'.format(NAME),
              sources = [WRAPPER_CPP] + CPP + UMICH_FILES,
              extra_compile_args = gccargs,
              libraries = libs,
              library_dirs = [],
              include_dirs = fixPath(["general", "klib", "vcf", "clusters", "libsrc", "merlin", "pdf"]) + ["source"]
              )
]
setup(name = NAME,
    version = VERSION,
    description = "A novel approach to use sequence data for linkage analysis",
    author = "Gao Wang",
    packages = [NAME],
    scripts = ['source/slinco'],
    package_dir = {NAME:'source'},
    cmdclass = {'build_py': build_py },
    ext_modules = LIBCORE_MODULE
)
