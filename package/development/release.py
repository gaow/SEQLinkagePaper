#!/usr/bin/env python2.7
# $File: release.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import os
import time
import sys
import subprocess
import shutil
import argparse
import zipfile
import platform

if sys.version.startswith('2.7.4') or sys.version.startswith('3.3.1'):
    sys.exit('Python 2.7.4 and 3.3.1 cannot be used due to a regression bug in the gzip module.')

###
SRCDIR = 'source'
EXEDIR = 'source'
DEVDIR = 'development'
EXE = ['seqlink']
PROJ = 'SEQLinkage'
VERDEV = '3.14159'
DEPS = ['tornado', 'brewer2mpl']
###

import imp
for item in DEPS:
    try:
        imp.find_module(item)
    except ImportError:
        sys.exit('Cannot build binary release: missing module "{}"!'.format(item))

def rm(x):
    # remove existing files or directories
    if os.path.isfile(x):
        os.remove(x)
    elif os.path.isdir(x):
        shutil.rmtree(x)

def modifyInit(version, include = [], exclude = []):
    # modify source/__init__.py to write version string
    try:
        rev = subprocess.check_output('date', shell=True).strip()
    except:
        rev = subprocess.check_output('date', shell=True).decode().strip()
    try:
        cppversion = subprocess.check_output('gcc -dumpversion', shell=True,
                                          stderr=subprocess.STDOUT).strip().split('\n')[0] 
    except:
        cppversion = subprocess.check_output('gcc -dumpversion', shell=True,
                                          stderr=subprocess.STDOUT).decode().strip().split('\n')[0] 
    full_version = "{}, released on {}, developed under Python {}.{}.{} and GNU/GCC {}".\
      format(version, rev, sys.version_info.major, sys.version_info.minor, sys.version_info.micro, cppversion)
    #
    content = []
    with open('{}/__init__.py'.format(SRCDIR), 'r') as init_file:
        for x in init_file.readlines():
            if x.startswith('VERSION'):
                content.append("VERSION = '{}'".format(version))
            elif x.startswith("FULL_VERSION"):
                content.append("FULL_VERSION = '{}'".format(full_version))
            else:
                if x.rstrip() not in include + exclude:
                    content.append(x.rstrip())
    content.extend(include)
    with open('{}/__init__.py'.format(SRCDIR), 'w') as init_file:
        init_file.write('\n'.join(content))

def setupEnvironment(version):
    for item in ['build', 'dist', '{}-{}'.format(PROJ, version)]:
        print('Removing directory {}'.format(item))
        rm(item)

def generateSWIGWrappers():
    # generate wrapper files to make sure they are up to date
    pass

def buildPackage(extra_args):
    try:
        print('Building and installing {} ...'.format(PROJ))
        with open(os.devnull, 'w') as fnull:
            ret = subprocess.call('python2.7 setup.py install ' + ' '.join(extra_args),
                                  shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to build and install {}.'.format(PROJ))
    except Exception as e:
        sys.exit('Failed to build and install {}: {}'.format(PROJ, e))


def buildSourcePackage(version):
    try:
        print('Building source package of {} {} ...'.format(PROJ, version))
        with open(os.devnull, 'w') as fnull:
            ret = subprocess.call('python2.7 setup.py sdist', shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to build source package of {}.'.format(PROJ))
        os.rename('dist/{}-{}.tar.gz'.format(PROJ, version),
            'dist/{}-{}-src.tar.gz'.format(PROJ, version).lower())
    except Exception as e:
        sys.exit('Failed to build source package of {}: {}'.format(PROJ, e))

def obtainPyInstaller(pyinstaller_dir):
    # check if pyinstaller is available
    #   if not, use git clone to get it
    #   if yes, try to update to the newest version
    pyinstaller_dir = os.path.expanduser(pyinstaller_dir.rstrip('/'))
    if pyinstaller_dir.endswith('pyinstaller'): pyinstaller_dir = pyinstaller_dir[:-12]
    if not os.path.isdir(pyinstaller_dir): pyinstaller_dir = os.getcwd()
    git_dir = os.path.join(pyinstaller_dir, 'pyinstaller')
    curdir = os.getcwd()
    if not os.path.isdir(git_dir):
        try:
            print('Downloading pyinstaller...')
            with open(os.devnull, 'w') as fnull:
                ret = subprocess.call('git clone git://github.com/pyinstaller/pyinstaller.git {}'.format(git_dir), shell=True, stdout=fnull)
                if ret != 0:
                    sys.exit('Failed to clone pyinstaller. Please check if you have git installed.'
                        'You can also get pyinstaller manually anf decompress it under the pyinstaller directory.')
        except Exception as e:
            sys.exit('Failed to clone pyinstaller: {}'.format(e) + 
                'You can get pyinstaller manually anf decompress it under the pyinstaller directory '
                'if you are having trouble getting git installed.')
    else:
        os.chdir(git_dir)
        try:
            print('Updating pyinstaller ...')
            with open(os.devnull, 'w') as fnull:
                ret = subprocess.call('git pull', shell=True, stdout=fnull)
                if ret != 0:
                    print('Failed to get latest version of pyinstaller. Using existing version.')
        except Exception as e:
            print('Failed to get latest version of pyinstaller ({}). Using existing version.'.format(e))
    os.chdir(curdir)
    return git_dir

def buildExecutables(git_dir, option):
    # use pyinstaller to create executable
    for exe in EXE:
        try:
            rm(os.path.join('dist', exe))
            print('Building executable {} ...'.format(exe))
            with open(os.devnull, 'w') as fnull:
                ret = subprocess.call('python2.7 {} {} --log-level=ERROR {} '
                    .format(os.path.join(git_dir, 'pyinstaller.py'), option, os.path.join(EXEDIR, exe)),
                    shell=True, stdout=fnull)
                if ret != 0:
                    sys.exit('Failed to create executable for command {}'.format(exe))
        except Exception as e:
            sys.exit('Failed to create executable for command {}: {}'.format(exe, e))


def createPackage(version, data_dir = None):
    machine = "MacOSX" if platform.system() == 'Darwin' else platform.system()
    bundle = '{}-{}-{}-{}'.format(PROJ, version, machine, platform.machine()).lower()
    dest = os.path.join('dist', bundle)
    rm(dest)
    os.makedirs(dest)
    # merge the installations
    for folder in ['dist/{}'.format(exe) for exe in EXE]:
        os.system('rsync -auqc {0}/* {1}'.format(folder, dest))
        rm(folder)
    # add data
    if data_dir:
        os.system('rsync -auqc {0}/* {1}'.format(data_dir, dest))
    with open('{}/installer/extractor-template.sh'.format(DEVDIR), 'r') as f:
        content = ["BUNDLE={0}\n".format(bundle) if x.startswith("BUNDLE=") else x for x in f.readlines()]
    with open(dest + '.bundle', 'w') as f:
        f.write(''.join(content))
    shutil.copy('{}/installer/install.sh'.format(DEVDIR), dest)
    os.system('tar c{1}f - -C dist {0} >> dist/{0}.bundle'.\
              format(bundle, 'j' if platform.system() == 'Darwin' else 'J'))
    rm(dest)

def createZipPackage(version):
    # after the creation of commands, create a zip file with OS and version information
    machine = "MacOSX" if platform.system() == 'Darwin' else platform.system()
    zipfilename = os.path.join('dist', '{}-{}-{}-{}.zip'
        .format(PROJ, version, machine, platform.machine()).lower())
    print('Adding executable to file {}'.format(zipfilename))
    with zipfile.ZipFile(zipfilename, 'w') as dist_file:
        for exe in EXE:
            cmd = '{}.exe'.format(exe) if os.name == 'win32' else exe
            dist_file.write(os.path.join('dist', cmd), cmd)
            rm(os.path.join('dist', cmd))

def tagRelease(version):
    try:
        ret = subprocess.check_output('git diff', shell=True)
        if ret:
            print('Commit all changes for the release of {}'.format(version))
            subprocess.call('git commit -am "Commit all change for the release of {}"'
                .format(version), shell=True)
        with open(os.devnull, 'w') as fnull:
            print('Tagging release {}...'.format(version))
            ret = subprocess.call('git tag -a {} -m "Version {} released at {}"; git push --tags'
                                  .format(version, version, time.asctime()), shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to tag release {}.'.format(version))
    except Exception as e:
        sys.exit('Failed to commit and tag release {}: {}'.format(version, e))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Create source distribution 
        and executables for a {} release. In addition to optional
        parameters version and tag, extra parameters would be specified and 
        will be passed directly to the 'python2.7 setup.py install' process.'''.format(PROJ))
    parser.add_argument('--version', default = VERDEV,
        help='''Modify {}/__init__.py to the specified version string and
            make the release.'''.format(SRCDIR))
    parser.add_argument('--tag', action='store_true',
        help='If specified, tag this release.')
    parser.add_argument('-p', '--pyinstaller_dir', default = '.',
        help='Path to the directory where pyinstaller git clone is located.')
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument('-o', '--onefile', action='store_true',
    #     help='Build one file executable instead of installer')
    parser.add_argument('-r', '--rebuild', action='store_true',
        help=argparse.SUPPRESS)
    #
    # allow recognied parameters to be set to the build process
    args, argv = parser.parse_known_args()
    args.onefile = False
    #
    modifyInit(args.version)
    if args.rebuild:
        setupEnvironment(args.version)
        generateSWIGWrappers()
    buildPackage(argv)
    buildSourcePackage(args.version)
    git_dir = obtainPyInstaller(args.pyinstaller_dir)
    if args.onefile:
        buildExecutables(git_dir, '--onefile')
        createZipPackage(args.version)  
    else:
        buildExecutables(git_dir, '--onedir')
        createPackage(args.version, data_dir = os.path.join(DEVDIR, 'data'))
    # if everything is OK, tag the release
    if args.tag:
        tagRelease(args.version)
    # if everything is done
    print('Source packages and executables are successfully generated and '
        'saved to directory dist.')
