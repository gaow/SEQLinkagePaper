# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import sys, imp, os
for item in ['faulthandler', 'numpy', 'matplotlib', 'prettyplotlib']:
    try:
        imp.find_module(item)
    except ImportError:
        sys.exit('Cannot build package: missing module "{}"!'.format(item))

try:
    imp.find_module('cstatgen')
except ImportError:
    from source.Utils import downloadURL
    from distutils.dir_util import remove_tree
    import tempfile
    cstatgen_url = "http://tigerwang.org/uploads/cstatgen.tar.gz"
    download_dir = tempfile.gettempdir()
    pkg = os.path.join(download_dir, "cstatgen.tar.gz")
    pkg_dir = os.path.join(download_dir, "cstatgen")
    sys.stderr.write("Downloading cstatgen library ...\n")
    downloadURL(cstatgen_url, download_dir, force = True)
    sys.stderr.write("Installing cstatgen library ...\n")
    try:
        remove_tree(pkg_dir)
    except:
        pass
    os.system("tar zxvf {} -C {} > /dev/null".format(pkg, download_dir))
    cwd = os.getcwd()
    os.chdir(pkg_dir)
    cmd = "python setup.py install {}".format(" ".join(sys.argv[2:]))
    os.system("{} > /dev/null".format(cmd))
    os.chdir(cwd)

from distutils.core import setup
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

from source import VERSION

NAME = "SEQLinco"

setup(name = NAME,
    version = VERSION,
    description = "A novel approach to use sequence data for linkage analysis",
    author = "Gao Wang",
    packages = [NAME],
    scripts = ['source/slinco'],
    package_dir = {NAME:'source'},
    cmdclass = {'build_py': build_py }
)
