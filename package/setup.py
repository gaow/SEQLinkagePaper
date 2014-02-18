# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import sys, imp
try:
    imp.find_module('cstatgen')
except ImportError:
    sys.exit('''Please install Gao Wang's "cstatgen" library''')

for item in ['faulthandler', 'numpy', 'matplotlib', 'prettyplotlib']:
    try:
        imp.find_module(item)
    except ImportError:
        sys.exit('Cannot build package: missing module "{}"!'.format(item))

from distutils.core import setup
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

from source import VERSION

NAME = "SEQLinco"

if sys.platform == "win32":
    sys.exit('Windows OS is not supported.')
    
setup(name = NAME,
    version = VERSION,
    description = "A novel approach to use sequence data for linkage analysis",
    author = "Gao Wang",
    packages = [NAME],
    scripts = ['source/slinco'],
    package_dir = {NAME:'source'},
    cmdclass = {'build_py': build_py }
)
