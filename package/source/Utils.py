import sys, os, subprocess, shutil, glob, shlex, urlparse, re, hashlib, tarfile
from multiprocessing import Process, Queue, Lock, Value
from distutils.dir_util import mkpath, remove_tree
from argparse import ArgumentParser, SUPPRESS
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

VERSION = "1.0.alpha"

class Environment:
    '''Ugly globals'''
    def __init__(self):
        self.__width_cache = 1
        # About the program 
        self.prog = 'slinco'
        self.proj = "SEQLinco"
        self.version = VERSION 
        # Runtime support
        self.resource_dir = os.path.expanduser('~/.{}'.format(self.proj))
        self.resource_bin = os.path.join(self.resource_dir, 'bin')
        self.cache_dir = os.path.join(os.getcwd(), 'cache')
        self.path = {'PATH':self.resource_bin}
        # File contents 
        self.build = 'hg19'
        self.delimiter = ' '
        self.ped_missing = ['0', '-9'] + ['none', 'null', 'na', 'nan', '.']
        self.missing = 'NA'
        self.trait = 'binary'
        # Input & output options
        self.output = 'LINKAGE'
        self.formats = {
            'plink':['.ped','.map'],
            'mega2':['.pre', '.map', '.name']
            }
        # Multiprocessing parameters
        self.batch = 50
        self.total_counter = Value('i',0)
        self.success_counter = Value('i',0)
        self.mendelerror_counter = Value('i',0)
            
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

###
# Utility functions
###

def flatten(l):
    return [item for sublist in l for item in sublist]

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

def downloadURL(URL, dest_dir, quiet = True, mode = None):
    if not os.path.isdir(dest_dir):
        mkpath(dest_dir)
    filename = os.path.split(urlparse.urlsplit(URL).path)[-1]
    dest = os.path.join(dest_dir, filename)
    if os.path.isfile(dest):
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

###
# Workhorse
###

def indexVCF(vcf):
    if not os.path.isfile(vcf + '.tbi'):
        env.log("Generating index file for [{}] ...".format(vcf))
        runCommand('tabix -p vcf {}'.format(vcf))
    return True

# @profile
def parseVCFline(line, exclude = [], info = None, letters = False, phased = True):
    if len(line) == 0:
        return []
    line = line.split('\t')
    v = ([line[3]] + line[4].split(',')) if letters else None
    gs = []
    for idx in range(len(line)):
        if idx < 9 or idx in exclude:
            continue
        if v:
            try:
                i, j = re.split('\/|\|', line[idx].split(":")[0])
                g = '{}{}'.format(v[int(i)], v[int(j)])
            except:
                try:
                    i = re.split('\/|\|', line[idx].split(":")[0])[0]
                    g = '{}'.format(v[int(i)])
                except:
                    g = '.'
        else:
            # Remove separater
            g = re.sub('\/|\|','',line[idx].split(":")[0])
        if not phased:
            g = ''.join(sorted(g))
        gs.append(g)
    return ([line[x-1] for x in info] + gs) if info is not None else gs
    
def extractSamplenames(vcf):
    samples = runCommand('tabix -H {}'.format(vcf)).strip().split('\n')[-1].split('\t')[9:]
    if not samples:
        env.error("Fail to extract samples from [{}]".format(vcf), exit = True)
    return samples

def extractFamilies(tfam):
    fam = {}
    samples = []
    with open(tfam, 'r') as f:
        for line in f.readlines():
            line = line.split()
            samples.append(line[1])
            if line[2] in env.ped_missing or line[3] in env.ped_missing:
                continue
            if line[0] not in fam:
                fam[line[0]] = [(line[2], line[3]), [line[1]]]
            else:
                fam[line[0]][1].append(line[1])
    for k in fam.keys():
        # remove missing father or mather
        if fam[k][0][0] not in samples or fam[k][0][1] not in samples:
            del fam[k]
            continue
        # remove parents from offsprings, just in case
        fam[k][1] = [x for x in fam[k][1] if x not in fam[k][0]]
    return fam

def rewriteFamfile(tfam, samples):
    '''remove samples from tfam file that are not in the
    input sample list'''
    if sorted(getColumn(tfam,2)) == sorted(samples):
        return
    os.rename(tfam, tfam + '.bak')
    try:
        with open(tfam + '.bak', 'r') as f:
            data = f.readlines()
        with open(tfam, 'w') as f:
            for item in data:
                if item.split()[1] in samples:
                    f.write(item)
    except:
        os.rename(tfam + '.bak', tfam)
        raise

def checkSamples(samp1, samp2):
    '''check if two sample lists agree
    1. samples in TFAM but not in VCF --> ERROR
    2. samples in VCF but not in TFAM --> give a message'''
    a_not_b = list(set(samp1).difference(set(samp2)))
    b_not_a = list(set(samp2).difference(set(samp1)))
    if b_not_a:
        raise RuntimeError('{:,d} samples found in TFAM file but not in VCF file:\n{}'.\
                           format(len(b_not_a), '\n'.join(b_not_a))) 
    if a_not_b:
        env.log('{:,d} samples in VCF file will be ignored due to absence in TFAM file'.format(len(a_not_b)))
    return True

def loadMap(maps):
    return {}

class Cache:
    def __init__(self, cache_dir, cache_name):
        self.cache_dir = cache_dir
        self.cache_name = os.path.join(cache_dir, cache_name + '.cache')
        self.cache_info = cache_info = os.path.join(cache_dir, '.info.' + cache_name)
        mkpath(cache_dir)

    def check(self):
        if not os.path.isfile(self.cache_info):
            return False
        with open(self.cache_info, 'r') as f:
            for item in f.readlines():
                name, key = item.split()
                if not os.path.isfile(name) or key != calculateFileMD5(name):
                    return False
        return True

    def load(self):
        with tarfile.open(self.cache_name) as f:
            f.extractall(self.cache_dir)        

    def write(self, ext = None, pre = None, otherfiles = []):
        '''Add files to cache'''
        if not self.check():
            with tarfile.open(self.cache_name, 'w:bz2') as f:
                for item in os.listdir(self.cache_dir):
                    if (item.endswith(ext) or ext is None) and (item.startswith(pre) or pre is None):
                        f.add(os.path.join(self.cache_dir, item), arcname=item)
            otherfiles.append(self.cache_name)
            signatures = ['{}\t{}'.format(x, calculateFileMD5(x)) for x in otherfiles]
            with open(self.cache_info, 'w') as f:
                f.write('\n'.join(signatures))

###
# Encoder core
###

class PseudoAutoRegion:
    def __init__(self, chrom, build):
        if build == ['hg18', 'build36'] and chrom.lower() in ['x', '23']:
            self.check = self.checkChrX_hg18
        elif build in ['hg18', 'build36'] and chrom.lower() in ['y', '24']:
            self.check = self.checkChrY_hg18
        elif build in ['hg19', 'build37'] and chrom.lower() in ['x', '23']:
            self.check = self.checkChrX_hg19
        elif build in ['hg19', 'build37'] and chrom.lower() in ['y', '24']:
            self.check = self.checkChrY_hg19
        else:
            self.check = self.notWithinRegion

    def checkChrX_hg18(self, pos):
        return (pos >= 1 and pos <= 2709520) or \
            (pos >= 154584238 and pos <= 154913754)

    def checkChrY_hg18(self, pos):
        return (pos >= 1 and pos <= 2709520) or \
            (pos >= 57443438 and pos <= 57772954)

    def checkChrX_hg19(self, pos):
        return (pos >= 60001 and pos <= 2699520) or \
            (pos >= 154931044 and pos <= 155270560)

    def checkChrY_hg19(self, pos):
        return (pos >= 10001 and pos <= 2649520) or \
            (pos >= 59034050 and pos <= 59373566)

    def notWithinRegion(self, pos):
        return False

class RData(dict):
    def __init__(self, samples, families):
        self.__families = families
        self.__families_flat = {key : flatten(families[key]) for key in families}
        self.__samples = samples
        self.reset()

    def reset(self):
        data = {}
        for item in self.__samples:
            data[item] = []
        data['__variant_info'] = []
        self.update(data)

    def samples(self):
        return self.__samples

    def families(self):
        # {fid : [(pid, mid), [sids ...]], ...}
        return self.__families

    def ffamilies(self):
        # {fid : [pid, mid, sids ...], ...}
        return self.__families_flat
    
class RegionExtractor:
    '''Extract given genomic region from VCF
    converting genotypes into dictionary of
    genotype list'''
    def __init__(self, filename, exclude_idx, build = env.build):
        self.vcf = filename
        self.exclude_idx = exclude_idx
        self.chrom = self.startpos = self.endpos = self.name = self.distance = None
        self.xchecker = PseudoAutoRegion('X', build)
        self.ychecker = PseudoAutoRegion('Y', build)
        
    def apply(self, data):
        # Only keep position information
        info = (2,)
        lines = [parseVCFline(x, exclude = self.exclude_idx, info=info) for x in \
                  runCommand('tabix {} {}:{}-{}'.\
                             format(self.vcf, self.chrom, self.startpos, self.endpos)).strip().split('\n')]
        if len(lines) == 0 or len(lines[0]) == 0:
            return False
        #
        if not len(lines[0]) - len(info) == len(data.samples()):
            raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name,len(lines[0]) - len(info), len(data.samples())))
        #
        for line in lines:
            data['__variant_info'].append(int(line[0]))
            # for each family, mark monomorphic site as missing
            # and assign their genotype 
            for k, eachfam in data.ffamilies().items():
                geno = [y for j, y in enumerate(line[1:]) if j in \
                        [i for i, x in enumerate(data.samples()) if x in eachfam]]
                if len(set(flatten([list(x) for x in set(geno) if x != '.']))) == 1:
                    # monomorphic site for this family
                    # skip this site
                    continue
                else:
                    for person, x in zip(eachfam, geno):
                        data[person].append(x)
                    data['__{}_major'.format(k)] = sorted(Counter(geno).most_common()[0][0])[0]
        return True
            
    def getRegion(self, region):
        self.chrom, self.startpos, self.endpos, self.name, self.distance = region
        if self.chrom.lower() in ['x','23','chrx','chr23']:
            if self.xchecker.check(int(self.startpos)) or self.xchecker.check(int(self.endpos)): 
                self.chrom = 'XY'
        if self.chrom.lower() in ['y','24','chry','chr24']:
            if self.ychecker.check(int(self.startpos)) or self.ychecker.check(int(self.endpos)):
                self.chrom = 'XY'        
    
class MendelianErrorChecker:
    '''Check and fix Mendelian error in
    offspring genotypes'''
    def __init__(self):
        self.lock = Lock()

    # @profile
    def apply(self, data):
        fam = data.families()
        for f in fam:
            map(lambda l: self.__check(data, fam[f][0], fam[f][1], l),
                range(len(data[fam[f][0][0]])))
        return True

    # @profile
    def __check(self, data, parents, kids, locus):
        """ check mendelian error for one marker
        input are parent ID's, offspring ids and locus index
        """
        gdad, gmom = data[parents[0]][locus], data[parents[1]][locus]
        if gdad == '.' and gmom == '.':
            return
        if gdad == '.':
            gdad = '01'
        if gmom == '.':
            gmom = '01'
        for k in kids:
            gkid = data[k][locus]
            if gkid == '.':
                # try to impute
                if gdad == '00' and gmom == '00':
                    data[k][locus] = '00'
                elif gdad == '11' and gmom == '11':
                    data[k][locus] = '11'       
                else:
                    continue
            else:
                try:
                    if (gkid[0] in gdad and gkid[1] in gmom) \
                      or (gkid[0] in gmom and gkid[1] in gdad):
                        continue
                    else:
                        with self.lock:
                            env.mendelerror_counter.value += 1
                        data[k][locus] = '.'
                except IndexError:
                    # kid is haploid, which means chrom X or Y variant on a boy
                    # FIXME: impossible to correctly identify unless chrom number is given
                    if gkid in gmom or gdad:
                        continue
                    else:
                        with self.lock:
                            env.mendelerror_counter.value += 1
                        data[k][locus] = '.'                        
        
class GenoEncoder:
    def __init__(self, wsize):
        self.size = wsize

    def apply(self, data):
        data['__variant_info'] = sum(data['__variant_info']) / len(data['__variant_info'])
        for k, eachfam in data.ffamilies().items():
            pool = []
            for person in eachfam:
                # set missing value to major allele
                if len(data[person]) == 0:
                    pool.extend(['0','0'])
                else:
                    pool.append(''.join([x[0] if x[0] != '.' else \
                                     data['__{}_major'.format(k)] for x in data[person]]))
                    pool.append(''.join([x[1] if len(x) > 1 else \
                                     (x[0] if x[0] != '.' else data['__{}_major'.format(k)]) \
                                     for x in data[person]]))
            # pool: [allele1@samp1, allele2@samp1, allele1@samp2, allele2@samp2, ...]
            if self.size > 1 or self.size == -1:
                pool = self.collapse(pool)
            # find unique pattern
            mapping = {item : idx + 1 for idx, item in enumerate(sorted(set(pool)))}
            # reset data
            idx = 0
            for person in eachfam:
                data[person] = '{} {}'.format(mapping[pool[idx]], mapping[pool[idx+1]])
                idx += 2
        # drop the region if it is monomorphic
        for k in data.samples():
            if data[k] != '1 1':
                # at least one sample is polymorphic
                return True
        return False
    
    def collapse(self, x):
        '''input data is a list of binary strings
        will be collapsed on every self.size characters'''
        # adjust size
        size = self.__adjust_size(len(x[0]))
        for idx, item in enumerate(x):
           x[idx] = ''.join(['1' if '1' in iitem else '0' \
                     for iitem in [item[i:i+size] for i in range(0, len(item), size)]])
        return x

    def __adjust_size(self, n):
        if self.size < 0:
            return n
        res, rem = divmod(n, self.size)
        # reduce size by x such that rem + res * x = self.size - x
        return self.size - (self.size - rem) / (res + 1)

class LinkageWritter:
    def __init__(self):
        self.lock = Lock()
        self.chrom = self.prev_chrom = self.name = self.distance = None
        self.__reset()

    def apply(self, data):
        # input data receieved, call it a "success"
        with self.lock:
            env.success_counter.value += 1
        if self.chrom != self.prev_chrom:
            if self.prev_chrom is None:
                self.prev_chrom = self.chrom
            else:
                # new chrom entered,
                # commit whatever is in buffer before accepting new data
                self.commit()
        # FIXME: currently only applies to collapsing theme
        self.output += env.delimiter.join(
            [self.chrom, self.name, self.distance, str(data['__variant_info'])] + \
                                          [data[s] for s in data.samples()]) + '\n'
        if self.counter < env.batch:
            self.counter += 1
        else:
            self.commit()

    def commit(self):
        if not self.output:
            return
        with self.lock:
            with open(os.path.join(env.cache_dir, '{}.chr{}.tped'.format(env.output, self.prev_chrom)),
                      'a') as f:
                f.write(self.output)
        self.__reset()

    def __reset(self):
        self.output = ''
        self.counter = 0
        self.prev_chrom = self.chrom
            
    def getRegion(self, region):
        self.chrom = region[0]
        self.name, self.distance = region[3:]

class EncoderWorker(Process):
    def __init__(self, wid, queue, data, rextractor, mchecker, gencoder, linkwritter):
        Process.__init__(self)
        self.worker_id = wid
        self.queue = queue
        self.modules = [rextractor, mchecker, gencoder, linkwritter]
        self.data = data
        self.lock = Lock()

    def report(self):
        env.log('{:,d} units processed; {:,d} Mendelian errors fixed; {:,d} super markers generated'.\
                format(env.total_counter.value, env.mendelerror_counter.value, env.success_counter.value),
                flush = True)
        
    def run(self):
        while True:
            try:
                region = self.queue.get()
                if region is None:
                    self.modules[-1].commit()
                    self.report()
                    break
                else:
                    with self.lock:
                        env.total_counter.value += 1
                    self.modules[0].getRegion(region)
                    self.modules[-1].getRegion(region)
                    self.data.reset()
                    for m in self.modules:
                        if not m.apply(self.data):
                            # previous module failed
                            break
                    if env.total_counter.value % (env.batch * env.jobs) == 0:
                        self.report()
            except KeyboardInterrupt:
                break

###
# Data converters
###

def formatPlink(tpeds, tfams, outdir):
    mkpath(outdir)
    cmds = ['transpose.pl --format plink --tped {} --tfam {} --out {}'.\
            format(i,j, os.path.join(outdir, os.path.splitext(os.path.split(i)[-1])[0])) for i, j in zip(tpeds, tfams)]
    env.jobs = max(min(args.jobs, cmds), 1)
    runCommands(cmds, env.jobs)
                
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
    env.jobs = max(min(args.jobs, cmds), 1)
    runCommands(cmds, env.jobs)


def formatMlinkfrommega2():
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

def runMlink():
    #if 'mlink' not in env.formats:
    #    formatMlink()
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    cmds = ['runMlink.pl MLINK/{}.{} {}'.format(env.output, chrs[i], env.resource_dir) for i in range(25)]
    env.jobs = max(min(args.jobs, cmds), 1)
    runCommands(cmds, env.jobs)

def plotMlink():
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    for dir in ['MLINK/{}.{}'.format(env.output, i) for i in chrs]:
        lods = []
        with open(dir + '/all_lodscores.txt', 'r') as f:
            for line in f.readlines():
                lod = line.split()[-1]
                lods.append(lod)
        lods = np.array(map(float,lods)).reshape((10,-1))
        plt.pcolormesh(lods)
        plt.savefig(dir + '/lods.pdf')
    #cmds = ['runMlink.pl MLINK/{}.{} {}'.format(env.output, chrs[i], env.resource_dir) for i in range(25)]
    #env.jobs = max(min(args.jobs, cmds), 1)
    #runCommands(cmds, env.jobs)
    
