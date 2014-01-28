# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

from distutils.core import setup, Extension
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py
import sys, os, subprocess
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
def fixPath(fn, prefix = "source/libmped"):
    if type(fn) is list:
        return [os.path.join(prefix, x) for x in fn]
    else:
        return os.path.join(prefix, fn)
#
WRAPPER_CPP = fixPath('chp_{0}.cpp'.format(PYVERSION))
WRAPPER_PY = fixPath('chp_{0}.py'.format(PYVERSION))
WRAPPER_I = fixPath('chp.i')
HEADER = fixPath(['chp.i', 'chp.hpp', 'Core.hpp', 'Exception.hpp'])
CPP = fixPath(['Core.cpp'])
# generate wrapper files
try:
    ret = subprocess.call(['swig', '-python', '-external-runtime', 'swigpyrun.h'], shell=False)
    if ret != 0:
        sys.exit('Failed to generate swig runtime header file. Please install swig.')
    #
    if (not os.path.isfile(WRAPPER_PY) or not os.path.isfile(WRAPPER_CPP) or \
        os.path.getmtime(WRAPPER_CPP) < max([os.path.getmtime(x) for x in HEADER + CPP])):
        ret = subprocess.call(['swig'] + SWIG_OPTS + ['-o', WRAPPER_CPP, WRAPPER_I], shell=False)
        if ret != 0:
            sys.exit('Failed to generate cpp extension.')
        os.rename('chp.py', WRAPPER_PY)
except OSError as e:
    sys.exit('Failed to generate wrapper file: {0}'.format(e))
#
MERLIN_FILES = glob(fixPath("clusters/*.cpp")) + \
  glob(fixPath("libsrc/*.cpp")) + \
  glob(fixPath("merlin/*.cpp")) + \
  glob(fixPath("pdf/*.cpp"))
# Under linux/gcc, lib stdc++ is needed for C++ based extension.
libs = ['stdc++'] if sys.platform == 'linux2' else []
gccargs = ["-O3", "-march=native", "-std=c++11"]
#
CHP_MODULE = [
    Extension('{}.libmped._chp'.format(NAME),
              sources = [WRAPPER_CPP] + CPP + MERLIN_FILES,
              extra_compile_args = gccargs,
              libraries = libs,
              library_dirs = [],
              include_dirs = fixPath([".", "clusters", "libsrc", "merlin", "pdf"])
              )
]
setup(name = NAME,
    version = VERSION,
    description = "A novel approach to use sequence data for linkage analysis",
    author = "Gao Wang",
    packages = [NAME, '{}.libmped'.format(NAME)],
    scripts = ['source/slinco'],
    package_dir = {NAME:'source'},
    cmdclass = {'build_py': build_py },
    ext_modules = CHP_MODULE
)
