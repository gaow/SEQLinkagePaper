#!/usr/bin/python
# Copyright (c) 2013 - 2014, Gao Wang <gaow@bcm.edu> and Di Zhang <di.zhang@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import sys, os, subprocess, shutil, glob, shlex, urlparse, re, hashlib, tempfile
from cStringIO import StringIO
from contextlib import contextmanager
from multiprocessing import Pool, Process, Queue, Lock, Value
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import prettyplotlib as ppl
from distutils.dir_util import mkpath, remove_tree
from . import VERSION

class Environment:
    '''Ugly globals'''
    def __init__(self):
        self.__width_cache = 1
        # About the program 
        self.proj = "SEQLinkage"
        self.prog = 'seqlink'
        self.version = VERSION 
        # Runtime support
        self.resource_dir = os.path.expanduser('~/.{}'.format(self.proj))
        self.resource_bin = os.path.join(self.resource_dir, 'bin')
        self.cache_dir = os.path.join(os.getcwd(), 'cache')
        self.tmp_dir = self.__mktmpdir()
        self.tmp_cache = os.path.join(self.tmp_dir, 'CACHE')
        self.path = {'PATH':"{}:{}".format(self.resource_bin, os.environ["PATH"])}
        self.debug = False
        self.quiet = False
        # File contents 
        self.build = 'hg19'
        self.delimiter = " "
        self.ped_missing = ['0', '-9'] + ['none', 'null', 'na', 'nan', '.']
        self.trait = 'binary'
        # Input & output options
        self.output = 'LINKAGE'
        self.tmp_log = os.path.join(self.tmp_dir, "clog." + self.output)
        # Multiprocessing counters
        self.batch = 50
        self.lock = Lock()
        self.total_counter = Value('i',0)
        self.success_counter = Value('i',0)
        self.null_counter = Value('i',0)
        self.trivial_counter = Value('i',0)
        self.chperror_counter = Value('i',0)
        self.variants_counter = Value('i',0)
        self.triallelic_counter = Value('i',0)
        self.commonvar_counter = Value('i',0)
        self.mendelerror_counter = Value('i',0)
        self.recomb_counter = Value('i',0)
        self.skipped_counter = Value('i',0)
        self.makeped_counter = Value('i',0)
        self.pedcheck_counter = Value('i',0)
        self.unknown_counter = Value('i',0)
        self.mlink_counter = Value('i',0)

    def __mktmpdir(self, where = None):
        class LockedTempDir(str):
            def __init__(self, path):
                self = path
                open(os.path.join(self, '.lock'), 'a').close()

            def __del__(self):
                try:
                    os.remove(os.path.join(self, '.lock'))
                except:
                    pass

        if where in [None, 'None', '']:
            where = tempfile.gettempdir()
        else:
            where = os.path.expanduser(where)
        if os.path.isdir(where) and ((not os.access(where, os.R_OK)) or (not os.access(where, os.W_OK))):
            self.error('Cannot set temporary directory to directory {} because '.format(where) + \
                       'it is not readable or writable.', exit = True)
        pattern = re.compile(r'{}_tmp_*(.*)'.format(self.proj))
        for fn in os.listdir(where):
            if pattern.match(fn) and not os.path.isfile(os.path.join(where, fn, '.lock')):
                try:
                    remove_tree(os.path.join(where, fn))
                except:
                    pass
        tmp = LockedTempDir(tempfile.mkdtemp(prefix='{}_tmp_'.format(self.proj), dir = where))
        mkpath(os.path.join(tmp, 'CACHE'))
        return tmp

    def ResetTempdir(self, path = None):
        self.tmp_dir = self.__mktmpdir(path)
        self.tmp_cache = os.path.join(self.tmp_dir, 'CACHE')
        self.tmp_log = os.path.join(self.tmp_dir, "clog." + self.output)
            
    def error(self, msg = None, show_help = False, exit = False):
        if msg is None:
            sys.stderr.write('\n')
            return
        if type(msg) is list:
            msg = ' '.join(map(str, msg))
        else:
            msg = str(msg)
        start = '\n' if msg.startswith('\n') else ''
        end = '\n' if msg.endswith('\n') else ''
        msg = msg.strip()
        sys.stderr.write(start + "\033[1;40;33mERROR: {}\033[0m\n".format(msg) + end)
        if show_help:
            self.log("Type '{} -h' for help message".format(env.prog))
            sys.exit()
        if exit:
            sys.exit()
        
    def log(self, msg = None, flush=False):
        if self.debug or self.quiet:
            return
        if msg is None:
            sys.stderr.write('\n')
            return
        if type(msg) is list:
            msg = ' '.join(map(str, msg))
        else:
            msg = str(msg)
        start = "{0:{width}}".format('\r', width = self.__width_cache + 10) + "\r" if flush else ''
        end = '' if flush else '\n'
        start = '\n' + start if msg.startswith('\n') else start
        end = end + '\n' if msg.endswith('\n') else end
        msg = msg.strip()
        if flush:
            self.__width_cache = len(msg)
        sys.stderr.write(start + "\033[1;40;32mMESSAGE: {}\033[0m".format(msg) + end)

env = Environment()

class StdoutCapturer(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

@contextmanager
def stdoutRedirect(to=os.devnull):
    '''
    import os

    with stdoutRedirect(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    fd = sys.stdout.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        sys.stdout.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w') # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'a') as file:
            _redirect_stdout(to=file)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

def runCommand(cmd, instream = None, msg = '', upon_succ = None, show_stderr = False, return_zero = True):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    popen_env = os.environ.copy()
    popen_env.update(env.path)
    try:
        tc = subprocess.Popen(cmd, stdin = subprocess.PIPE,
                              stdout = subprocess.PIPE, stderr = subprocess.PIPE,
                              env=popen_env)
        if instream:
            if sys.version_info.major == 3:
                instream = instream.encode(sys.getdefaultencoding())
            out, error = tc.communicate(instream)
        else:
            out, error = tc.communicate()
        if sys.version_info.major == 3:
            out = out.decode(sys.getdefaultencoding())
            error = error.decode(sys.getdefaultencoding())
        if return_zero:
            if tc.returncode < 0:
                raise ValueError ("Command '{0}' was terminated by signal {1}".format(cmd, -tc.returncode))
            elif tc.returncode > 0:
                raise ValueError ("{0}".format(error))
        if error and show_stderr:
            env.log(error)
    except OSError as e:
        raise OSError ("Execution of command '{0}' failed: {1}".format(cmd, e))
    # everything is OK
    if upon_succ:
        # call the function (upon_succ) using others as parameters.
        upon_succ[0](*(upon_succ[1:]))
    return out, error

class CMDWorker(Process):
    def __init__(self, queue):
        Process.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            try:
                cmd = self.queue.get()
                if cmd is None:
                    break
                else:
                    runCommand(cmd)
            except KeyboardInterrupt:
                break
            
def runCommands(cmds, ncpu):
    try:
        jobs = []
        queue = Queue()
        for i in cmds:
            queue.put(i)
        for i in range(ncpu):
            p = CMDWorker(queue)
            p.start()
            jobs.append(p)
            queue.put(None)
        for j in jobs:
            j.join()
    except KeyboardInterrupt:
        raise ValueError('Commands terminated!')

def downloadURL(URL, dest_dir, quiet = True, mode = None, force = False):
    if not os.path.isdir(dest_dir):
        mkpath(dest_dir)
    filename = os.path.split(urlparse.urlsplit(URL).path)[-1]
    dest = os.path.join(dest_dir, filename)
    if os.path.isfile(dest):
        if force:
            os.remove(dest)
        else:
            return
    # use wget
    try:
        # for some strange reason, passing wget without shell=True can fail silently.
        p = subprocess.Popen('wget {} -O {} {}'.format('-q' if quiet else '', dest, URL), shell=True)
        ret = p.wait()
        if ret == 0 and os.path.isfile(dest):
            if mode is not None:
                subprocess.Popen('chmod {} {}'.format(mode, dest), shell=True)
            return dest
        else:
            try:
                os.remove(dest)
            except OSError:
                pass
            raise RuntimeError('Failed to download {} using wget'.format(URL))
    except (RuntimeError, ValueError, OSError):
        # no wget command
        env.error('Failed to download {}'.format(filename))

def calculateFileMD5(filename):
    md5 = hashlib.md5()
    # limit the calculation to the first 1G of the file content
    block_size = 2**20  # buffer of 1M
    filesize = os.path.getsize(filename)
    try:
        if filesize < 2**26:
            # for file less than 1G, use all its content
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    if not data:
                        break
                    md5.update(data)
        else:
            count = 64
            # otherwise, use the first and last 500M
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    count -= 1
                    if count == 32:
                        f.seek(-2**25, 2)
                    if not data or count == 0:
                        break
                    md5.update(data)
    except IOError as e:
        sys.exit('Failed to read {}: {}'.format(filename, e))
    return md5.hexdigest()
 
def zipdir(path, zipfile, arcroot = '/'):
    path = os.path.normpath(path)
    for root, dirs, files in os.walk(path):
        for f in files:
            zipfile.write(os.path.join(root, f), arcname = os.path.join(arcroot, root[len(path) + 1:], f))

def removeFiles(dest, exclude = [], hidden = False):
    if os.path.isdir(dest):
        for item in os.listdir(dest):
            if item.startswith('.') and hidden == False:
                continue
            if os.path.splitext(item)[1] not in exclude:
                try:
                    os.remove(os.path.join(dest,item))
                except OSError:
                    pass

def copyFiles(pattern, dist, ignore_hidden = True):
    mkpath(dist)
    for fl in glob.glob(pattern):
        if os.path.isfile(fl):
            shutil.copy(fl, dist)

def downloadResources(fromto):
    for idx, item in enumerate(fromto):
        env.log('Checking local resources {0}/{1} ...'.format(idx + 1, len(fromto)), flush = True)
        downloadURL(item[0], item[1], mode = 777)
    env.log() 
    return True

def getColumn(fn, num, delim = None, exclude = None):
    num = num - 1
    with open(fn) as inf:
        output = []
        for line in inf:
            parts = line.split(delim) if delim is not None else line.split()
            if len(parts) > num and parts[num] != exclude:
                output.append(parts[num])
    return output

def wordCount(filename):
    """Returns a word/count dict for this filename."""
    word_count = {}
    with open(filename, 'r') as input_file:
        for line in input_file:
           words = line.split()
           for word in words:
               word = word.lower()
               if not word in word_count:
                   word_count[word] = 1
               else:
                   word_count[word] += 1
    return word_count

def connected_components(lists):
    neighbors = defaultdict(set)
    seen = set()
    for each in lists:
        for item in each:
            neighbors[item].update(each)
    def component(node, neighbors=neighbors, seen=seen, see=seen.add):
        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            see(node)
            nodes |= neighbors[node] - seen
            yield node
    for node in neighbors:
        if node not in seen:
            yield sorted(component(node))

def listit(t):
    return list(map(listit, t)) if isinstance(t, (list, tuple)) else t
            
def parseVCFline(line, exclude = []):
    if len(line) == 0:
        return None
    line = line.split('\t')
    # Skip tri-allelic variant
    if "," in line[4]:
        with env.lock:
            env.triallelic_counter.value += 1
        return None
    gs = []
    for idx in range(len(line)):
        if idx < 9 or idx in exclude:
            continue
        else:
            # Remove separater
            g = re.sub('\/|\|','',line[idx].split(":")[0])
            if g == '.' or g == '..':
                gs.append("00")
            else:
                gs.append(g.replace('1','2').replace('0','1'))
    return (line[0], line[1], line[3], line[4], line[2]), gs


###
# Check parameter input
###

def checkParams(args):
    '''set default arguments or make warnings'''
    env.debug = args.debug
    env.quiet = args.quiet
    args.vcf = os.path.abspath(os.path.expanduser(args.vcf))
    args.tfam = os.path.abspath(os.path.expanduser(args.tfam))
    for item in [args.vcf, args.tfam]:
        if not os.path.exists(item):
            env.error("Cannot find file [{}]!".format(item), exit = True)
    if args.output:
        env.output = os.path.split(args.output)[-1]
        env.tmp_log = os.path.join(env.tmp_dir, env.output + ".STDOUT")
    #
    if len([x for x in set(getColumn(args.tfam, 6)) if x.lower() not in env.ped_missing]) > 2:
        env.trait = 'quantitative'
    env.log('{} trait detected in [{}]'.format(env.trait.capitalize(), args.tfam))
    if not args.blueprint:
        args.blueprint = os.path.join(env.resource_dir, 'genemap.txt')
    args.format = [x.lower() for x in set(args.format)]
    if args.run_linkage and "linkage" not in args.format:
        args.format.append('linkage')
    if None in [args.inherit_mode, args.prevalence, args.wild_pen, args.muta_pen] and "linkage" in args.format:
        env.error('To generate LINKAGE format or run LINKAGE analysis, please specify all options below:\n\t--prevalence, -K\n\t--moi\n\t--wild-pen, -W\n\t--muta-pen, -M', show_help = True, exit = True)
    if args.tempdir is not None:
        env.ResetTempdir(args.tempdir)
    return True

###
# Run External tools
###

def indexVCF(vcf, verbose = True):
    if not vcf.endswith(".gz"):
        if os.path.exists(vcf + ".gz"):
            if verbose:
                env.error("Cannot compress [{0}] because [{0}.gz] exists!".format(vcf), exit = True)
            else:
                sys.exit()
        if verbose:
            env.log("Compressing file [{0}] to [{0}.gz] ...".format(vcf))
        runCommand('bgzip {0}'.format(vcf))
        vcf += ".gz"
    if not os.path.isfile(vcf + '.tbi') or os.path.getmtime(vcf) > os.path.getmtime(vcf + '.tbi'):
        if verbose:
            env.log("Generating index file for [{}] ...".format(vcf))
        runCommand('tabix -p vcf -f {}'.format(vcf))
    return vcf
    
def extractSamplenames(vcf):
    samples = runCommand('tabix -H {}'.format(vcf))[0].strip().split('\n')[-1].split('\t')[9:]
    if not samples:
        env.error("Fail to extract samples from [{}]".format(vcf), exit = True)
    return samples
