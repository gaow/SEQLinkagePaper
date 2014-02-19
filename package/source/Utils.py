#!/usr/bin/python
# Copyright (c) 2013, Gao Wang <gaow@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

import sys, os, subprocess, shutil, glob, shlex, urlparse, re, hashlib, tarfile, tempfile
from cStringIO import StringIO
from contextlib import contextmanager
from multiprocessing import Pool, Process, Queue, Lock, Value
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
        self.proj = "SEQLinco"
        self.prog = 'slinco'
        self.version = VERSION 
        # Runtime support
        self.resource_dir = os.path.expanduser('~/.{}'.format(self.proj))
        self.resource_bin = os.path.join(self.resource_dir, 'bin')
        self.cache_dir = os.path.join(os.getcwd(), 'cache')
        self.tmp_dir = self.__mktmpdir()
        self.path = {'PATH':self.resource_bin}
        self.debug = False
        # File contents 
        self.build = 'hg19'
        self.delimiter = " "
        self.ped_missing = ['0', '-9'] + ['none', 'null', 'na', 'nan', '.']
        self.trait = 'binary'
        # Input & output options
        self.output = 'LINKAGE'
        self.outputfam = os.path.join(self.cache_dir, '{}.tfam'.format(self.output))
        self.tmp_log = os.path.join(self.tmp_dir, self.output)
        self.formats = {
            'plink':['.ped','.map'],
            'mega2':['.pre', '.map', '.name']
            }
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
        self.mendelerror_counter = Value('i',0)
        self.recomb_counter = Value('i',0)

    def __mktmpdir(self):
        pattern = re.compile(r'{}_tmp_*(.*)'.format(self.proj))
        for fn in os.listdir(tempfile.gettempdir()):
            if pattern.match(fn):
                remove_tree(os.path.join(tempfile.gettempdir(), fn))
        return tempfile.mkdtemp(prefix='{}_tmp_'.format(self.proj))
            
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
            remove_tree(self.tmp_dir)
            sys.exit()
        if exit:
            remove_tree(self.tmp_dir)
            sys.exit()
        
    def log(self, msg = None, flush=False):
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

def runCommand(cmd, instream = None, msg = '', upon_succ=None):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    try:
        tc = subprocess.Popen(cmd, stdin = subprocess.PIPE,
                              stdout = subprocess.PIPE, stderr = subprocess.PIPE,
                              env=env.path)
        if instream:
            if sys.version_info.major == 3:
                instream = instream.encode(sys.getdefaultencoding())
            out, error = tc.communicate(instream)
        else:
            out, error = tc.communicate()
        if sys.version_info.major == 3:
            out = out.decode(sys.getdefaultencoding())
            error = error.decode(sys.getdefaultencoding())
        if tc.returncode < 0:
            raise ValueError ("Command '{0}' was terminated by signal {1}".format(cmd, -tc.returncode))
        elif tc.returncode > 0:
            raise ValueError ("{0}".format(error))
        else:
            if error:
                env.log(error)
    except OSError as e:
        raise OSError ("Execution of command '{0}' failed: {1}".format(cmd, e))
    # everything is OK
    if upon_succ:
        # call the function (upon_succ) using others as parameters.
        upon_succ[0](*(upon_succ[1:]))
    return out

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
 
def removeFiles(dest, exclude = [], hidden = False):
    if os.path.isdir(dest):
        for item in os.listdir(dest):
            if item.startswith('.') and hidden == False:
                continue
            if os.path.splitext(item)[1] not in exclude:
                os.remove(os.path.join(dest,item))

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
    args.vcf = os.path.abspath(os.path.expanduser(args.vcf))
    args.tfam = os.path.abspath(os.path.expanduser(args.tfam))
    for item in [args.vcf, args.tfam]:
        if not os.path.exists(item):
            env.error("Cannot find file [{}]!".format(item), exit = True)
    if args.output:
        env.output = os.path.split(args.output)[-1]
        env.outputfam = os.path.join(env.cache_dir, '{}.tfam'.format(env.output))
        env.tmp_log = os.path.join(env.tmp_dir, env.output)
    #
    if len([x for x in set(getColumn(args.tfam, 6)) if x.lower() not in env.ped_missing]) > 2:
        env.trait = 'quantitative'
    env.log('{} trait detected in [{}]'.format(env.trait.capitalize(), args.tfam))
    if not args.blueprint:
        args.blueprint = os.path.join(env.resource_dir, 'genemap.txt')
    # pop plink/mega2 format to first
    args.format = [x.lower() for x in set(args.format)]
    for item in ['mega2', 'plink']:
        if item in args.format:
            args.format.insert(0, args.format.pop(args.format.index(item)))
    return True

###
# Run External tools
###

def indexVCF(vcf):
    if not vcf.endswith(".gz"):
        if os.path.exists(vcf + ".gz"):
            env.error("Cannot compress [{0}] because [{0}.gz] exists!".format(vcf), exit = True)
        env.log("Compressing file [{0}] to [{0}.gz] ...".format(vcf))
        runCommand('bgzip {0}'.format(vcf))
        vcf += ".gz"
    if not os.path.isfile(vcf + '.tbi') or os.path.getmtime(vcf) > os.path.getmtime(vcf + '.tbi'):
        env.log("Generating index file for [{}] ...".format(vcf))
        runCommand('tabix -p vcf -f {}'.format(vcf))
    return vcf
    
def extractSamplenames(vcf):
    samples = runCommand('tabix -H {}'.format(vcf)).strip().split('\n')[-1].split('\t')[9:]
    if not samples:
        env.error("Fail to extract samples from [{}]".format(vcf), exit = True)
    return samples


def formatPlink(tpeds, tfams, outdir):
    mkpath(outdir)
    cmds = ['transpose.pl --format plink --tped {} --tfam {} --out {}'.\
            format(i,j, os.path.join(outdir, os.path.splitext(os.path.split(i)[-1])[0])) for i, j in zip(tpeds, tfams)]
    runCommands(cmds, max(min(env.jobs, cmds), 1))
                
def formatMega2(plinkfiles):
    copyFiles('PLINK/{}*'.format(env.output), 'MEGA2')
    trait = 'A' if env.trait == 'binary' else 'T'
    preheader = ['Pedigree', 'ID', 'Father', 'Mother', 'Sex', 'Trait.{}'.format(trait)]
    maps = loadMap(os.path.join(env.resource_dir, 'genemap.pkl')) 
    for item in sorted(glob.glob(plinkfiles), key=lambda x: os.path.splitext(x)[1], reverse = True):
        if item.endswith('.ped'):
            base = os.path.splitext(item)[0]
            os.rename(item, base + '.pre')
            # markers = ['M{0}.{1}.1 M{0}.{1}.2'.format(i+1, j) for i, j in enumerate([x for x in getColumn(base + '.map', 1) if x != 'Chromosome'])]
            marker_names = getColumn(base + '.map', 2)
            markers = ['{0}.M.1 {0}.M.1'.format(m) for m in marker_names]
            with open(base + '.pre', 'r+') as f:
                data = f.read()
                f.seek(0)
                f.write(env.delimiter.join(preheader + markers) + '\n' + data)
            with open(base + '.name', 'w') as f:
                f.write(env.delimiter.join(['Type', 'Name']) + '\n')
                f.write(env.delimiter.join([trait, 'Trait']) + '\n')
                for m in marker_names:
                    f.write(env.delimiter.join(['M', m]) + '\n')
        elif item.endswith('.map'):
            with open(item) as f:
                data = f.readlines()
            with open(item + '.tmp', 'w') as f:
                f.write(env.delimiter.join(['Chromosome', 'Map.k.a', 'Name', 'Map.k.m', 'Map.k.f', '{}.p'.format(env.build)]) + '\n')
                for line in data:
                    line = line.split()
                    if line[1] in maps:
                        f.write(env.delimiter.join([line[0], maps[line[1]][0], line[1], maps[line[1]][1], maps[line[1]][2], line[3]]) + '\n')
                    else:
                        # use original distance 
                        f.write(env.delimiter.join([line[0], line[2], line[1], line[2], line[2], line[3]]) + '\n')
            os.rename(item + '.tmp', item)
        else:
           continue

def formatMlink(tpeds, tfams, outdir):
    mkpath(outdir)
    cmds = ['transpose.pl --format mlink --tped {} --tfam {} --out {}'.\
            format(i,j, os.path.join(outdir, os.path.splitext(os.path.split(i)[-1])[0])) for i, j in zip(tpeds, tfams)]
    runCommands(cmds, max(min(env.jobs, cmds), 1))


def formatMlinkformmega2():
    if 'mega2' not in env.formats:
	mkpath('MEGA2')
        formatMega2('MEGA2/{}*'.format(env.output))
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    #template = os.path.join(env.resource_dir, 'MEGA2.template')
    for chrNum in range(25):
        pedfile = 'MEGA2/{}.{}.pre'.format(env.output, chrs[chrNum])
        mapfile = 'MEGA2/{}.{}.map'.format(env.output, chrs[chrNum])
        namefile = 'MEGA2/{}.{}.name'.format(env.output, chrs[chrNum])
        if not (os.path.isfile(pedfile) and os.path.isfile(mapfile) and os.path.isfile(namefile)):
            continue
        outpath = 'MLINK/{}.{}'.format(env.output,chrs[chrNum])
        batchfile = os.path.join(outpath, 'MEGA2.BATCH')
        mkpath(outpath)
        with open(batchfile, 'w') as f:
            f.write('Input_Pedigree_File={}'.format(pedfile) + '\n')
            f.write('Input_Locus_File={}'.format(namefile) + '\n')
            f.write('Input_Map_File={}'.format(mapfile) + '\n')
            f.write('Input_Untyped_Ped_Option=2' + '\n')
            f.write('Input_Do_Error_Sim=no' + '\n')
            f.write('Output_Path={}'.format(outpath) + '\n')
            f.write('AlleleFreq_SquaredDev=999999999.000000' + '\n')
            #f.write('Value_Marker_Compression=1' + '\n')
            f.write('Analysis_Option=Vitesse' + '\n')
            f.write('Count_Genotypes=1' + '\n')
            f.write('Count_Halftyped=no' + '\n')
            f.write('Value_Genetic_Distance_Index=0' + '\n')
            f.write('Value_Genetic_Distance_SexTypeMap=0' + '\n')
            f.write('Value_Base_Pair_Position_Index=1' + '\n')
            f.write('Chromosome_Single={}'.format(chrNum + 1) + '\n')
            f.write('Traits_Combine=1 2 e' + '\n')
            f.write('Analysis_Sub_Option=MLINK' + '\n')
            f.write('Default_Outfile_Names=yes' + '\n')
        runCommand('mega2 -wx {}'.format(batchfile))

def runMlink(blueprint):
    #if 'mlink' not in env.formats:
    #    formatMlink()
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    cmds = ['runMlink.pl MLINK/{}.{} {} {}'.format(env.output, chrs[i], env.resource_dir, blueprint) for i in range(25)]
    runCommands(cmds, max(min(env.jobs, cmds), 1))

def heatmap(dir):
    env.log("start ploting heatmap for" + dir + "\n")
    lods = []
    with open(dir + '/all_lodscores.txt', 'r') as f:
        for line in f.readlines():
            lod = line.split()[-1]
            lods.append(lod)
        lods = np.array(map(float,lods)).reshape((10,-1))
        ppl.pcolormesh(lods)
        plt.savefig('MLINK/heatmap/{}.lods.png'.format(os.path.basename(dir)))
    env.log("end ploting heatmap for" + dir + "\n")

def plotMlink():
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    dirs = filter(lambda x: os.path.exists(x), ['MLINK/{}.{}'.format(env.output, i) for i in chrs])
    mkpath('MLINK/heatmap')
    pool = Pool(env.jobs)
    pool.map(heatmap, dirs)
