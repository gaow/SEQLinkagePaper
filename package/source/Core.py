#!/usr/bin/python
# Copyright (c) 2013, Gao Wang <gaow@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

from SEQLinco.Utils import *
from multiprocessing import Process, Queue
from collections import Counter, OrderedDict
from itertools import chain
import sys, faulthandler

if sys.version_info.major == 2:
    from cstatgen import cstatgen_py2 as cstatgen
else:
    from cstatgen import cstatgen_py3 as cstatgen

def checkVCFBundle(vcf):
    '''VCF bundle should have a .gz file and a tabix file'''
    if not vcf.endswith(".gz"):
        env.error("Input VCF file has to be bgzipped and indexed (http://samtools.sourceforge.net/tabix.shtml).",
                  exit = True)
    if not os.path.isfile(vcf + '.tbi'):
        env.error("Index file [{}] not found.".format(vcf + '.tbi'), exit = True)
    return True

def rewriteFamfile(tfam, samples, keys):
    with open(tfam, 'w') as f:
        f.write('\n'.join(['\t'.join(samples[k]) for k in keys]))
        f.write('\n')

def checkSamples(samp1, samp2):
    '''check if two sample lists agree
    1. samples in TFAM but not in VCF --> ERROR
    2. samples in VCF but not in TFAM --> give a message'''
    a_not_b = list(set(samp1).difference(set(samp2)))
    b_not_a = list(set(samp2).difference(set(samp1)))
    if b_not_a:
        env.error('{:,d} samples found in TFAM file but not in VCF file:\n{}'.\
                           format(len(b_not_a), '\n'.join(b_not_a)), exit = True) 
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

    def write(self, ext = None, pre = None, files = [], otherfiles = []):
        '''Add files to cache'''
        if not self.check():
            with tarfile.open(self.cache_name, 'w:bz2') as f:
                for item in os.listdir(self.cache_dir):
                    if ((item.endswith(ext) or ext is None) and (item.startswith(pre) or pre is None)) \
                      or item in files:
                        f.add(os.path.join(self.cache_dir, item), arcname=item)
            otherfiles.append(self.cache_name)
            signatures = ['{}\t{}'.format(x, calculateFileMD5(x)) for x in otherfiles]
            with open(self.cache_info, 'w') as f:
                f.write('\n'.join(signatures))

    def clear(self, pre, ext):
        for fl in glob.glob(os.path.join(self.cache_dir, pre) +  "*" + ext) + [self.cache_info, self.cache_name]:
            if os.path.isfile(fl):
                os.remove(fl)
        

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

class TFAMParser:
    def __init__(self, tfam):
        self.families, self.samples = self.parse(tfam)

    def __add_or_app(self, obj, key, value):
        islist = type(value) is list
        if key not in obj:
            if not islist:
                obj[key] = [value]
            else:
                obj[key] = value
        else:
            if value not in obj[key]:
                if not islist:
                    obj[key].append(value)
                else:
                    obj[key].extend(value)

    def parse(self, tfam):
        '''Rules:
        1. samples have to have unique names
        2. both parents for a sample should be available
        3. founders should have at least one offspring'''
        fams = {}
        observedFounders = {}
        expectedFounders = {}
        samples = {}
        with open(tfam, 'r') as f:
            for idx, line in enumerate(f.readlines()):
                line = line.split()
                if len(line) != 6:
                    env.error("skipped line {} (has {} != 6 columns!)".format(idx, len(line)))
                    continue
                if line[1] in samples:
                    env.error("skipped line {} (duplicate sample name '{}' found!)".format(idx, line[1]))
                    continue
                # collect sample info
                samples[line[1]] = [line[0], line[1], line[2], line[3], line[4], line[5]]
                # collect family member
                self.__add_or_app(fams, line[0], line[1])
                # collect founders for family
                if line[2] in env.ped_missing and line[3] in env.ped_missing:
                    self.__add_or_app(observedFounders, line[0], line[1])
                else:
                    self.__add_or_app(expectedFounders, line[0], (line[1], line[2], line[3]))
        # if a sample's expected founder is not observed in tfam
        # then the sample itself is a founder 
        # after this, all families should have founders
        for k in expectedFounders.keys():
            if k not in observedFounders:
                for item in expectedFounders[k]:
                    samples[item[0]][2] = samples[item[0]][3] = "0"
                observedFounders[k] = [item[0] for item in expectedFounders[k]]
                continue
            for item in expectedFounders[k]:
                if (item[1] not in observedFounders[k] and item[2] in observedFounders[k]) \
                    or (item[2] not in observedFounders[k] and item[1] in observedFounders[k]) \
                    or (item[1] not in observedFounders[k] and item[2] not in observedFounders[k]):
                    # missing one or two founders
                    samples[item[0]][2] = samples[item[0]][3] = "0"
                    observedFounders[k].append(item[0])
        # now remove trivial families 
        for k in fams.keys():
            if Counter(observedFounders[k]) == Counter(fams[k]):
                del fams[k]
                continue
        #
        valid_samples = []
        for value in fams.values():
            valid_samples.extend(value)
        samples = {k : samples[k] for k in valid_samples}
        return fams, samples
    
class RData(dict):
    def __init__(self, samples, families, samples_vcf):
        # a dict of {sid:[fid, pid, mid, sex, trait], ...}
        self.samples = OrderedDict()
        for k in samples_vcf:
            if k in samples:
                self.samples[k] = samples[k]
        # a dict of {fid:[member names], ...}
        self.families = {}
        # a dict of {fid:[idx ...], ...}
        self.famsampidx = {}
        # a dict of {fid:[idx ...], ...}
        self.famvaridx = {}
        # reorder family samples based on order of VCF file
        for k in families:
            self.families[k] = [x for x in self.samples if x in families[k]]
            self.famsampidx[k] = [i for i, x in enumerate(samples_vcf) if x in families[k]]
        self.reset()

    def reset(self):
        for item in self.samples:
            self[item] = []
        self.variants = []
        self.chrom = None
        for k in self.families:
            self.famvaridx[k] = []
        self.superMarkerCount = 0

    def getMidPosition(self):
        if len(self.variants) == 0:
            return None
        return sum([x[1] for x in self.variants]) / len(self.variants)

    def getFamVariants(self, fam, style = None):
        famvar = [item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam]]
        if style is None:
            return famvar
        elif style == "map":
            names = []
            pos = []
            for idx, item in enumerate(famvar):
                # names.append("{}-{}-{}".format(item[-1], item[1], idx))
                names.append("V{}".format(idx))
                pos.append(item[1])
            return names, pos
        else:
            raise ValueError("Unknown style '{}'".format(style))

    def getFamSamples(self, fam):
        output = [[]] * len(self.families[fam])
        for idx, item in enumerate(self.families[fam]):
            output[idx] = self.samples[item][:-1]
            output[idx].extend(self[item])
        # Input for CHP requires founders precede offsprings!
        return sorted(output, key = lambda x: x[2])


class RegionExtractor:
    '''Extract given genomic region from VCF
    converting genotypes into dictionary of
    genotype list'''
    def __init__(self, filename, build = env.build):
        self.vcf = cstatgen.VCFstream(filename)
        self.chrom = self.startpos = self.endpos = self.name = self.distance = None
        self.xchecker = PseudoAutoRegion('X', build)
        self.ychecker = PseudoAutoRegion('Y', build)
        
    def apply(self, data):
        # Clean up
        data.reset()
        data.chrom = self.chrom
        self.vcf.Extract(self.chrom, self.startpos, self.endpos)
        lineCount = 0
        # for each variant site
        while (self.vcf.Next()):
            # skip tri-allelic sites
            if not self.vcf.IsBiAllelic():
                with env.lock:
                    env.triallelic_counter.value += 1
                continue
            # check if the line's sample number matches the entire VCF sample number 
            if not self.vcf.CountSampleGenotypes() == self.vcf.sampleCount:
                raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name, self.vcf.CountSampleGenotypes(), self.vcf.sampleCount))
            # valid line found
            data.variants.append([self.vcf.GetChrom(), self.vcf.GetPosition(), self.name])
            # for each family assign member genotype if the site is non-trivial in the fam
            for k in data.families:
                gs = self.vcf.GetGenotypes(data.famsampidx[k])
                if len(set(''.join([x for x in gs if x != "00"]))) == 1:
                    # skip monomorphic gs
                    continue
                else:
                    # this variant is found in the family
                    data.famvaridx[k].append(lineCount)
                    for person, g in zip(data.families[k], gs):
                        data[person].append(g)
            lineCount += 1
        # 
        if lineCount == 0:
            return 1
        else:
            with env.lock:
                env.variants_counter.value += lineCount 
            return 0
            

    def getRegion(self, region):
        self.chrom, self.startpos, self.endpos, self.name, self.distance = region
        self.chrom = self.chrom.upper().replace("CHR", "")
        self.startpos = int(self.startpos)
        self.endpos = int(self.endpos)
        if self.chrom in ['X','23']:
            if self.xchecker.check(self.startpos) or self.xchecker.check(self.endpos): 
                self.chrom = 'XY'
        if self.chrom in ['Y','24']:
            if self.ychecker.check(self.startpos) or self.ychecker.check(self.endpos):
                self.chrom = 'XY'

    
class MarkerMaker:
    def __init__(self, wsize):
        self.size = wsize
        # adjust physical distance to map distance, 1 / 100 million
        self.position_adj = 1 / float(100000000)
        self.missings = ("0", "0")
        self.missing = "0"

    def apply(self, data):
        # data.superMarkerCount is the max num. of recombinant fragments among all fams
        for item in data.families:
            try:
                varnames, varpos = data.getFamVariants(item, style = "map")
                if len(varnames) < 2:
                    # no more than 2 variants found in family
                    for person in data.families[item]:
                        if len(data[person]) == 0:
                            data[person] = self.missings
                        else:
                            data[person] = tuple(data[person][0])
                    if data.superMarkerCount < 1:
                        data.superMarkerCount = 1
                else:
                    worker = cstatgen.CHP(self.size, self.position_adj, env.debug)
                    with stdoutRedirect(to = env.tmp_log + str(os.getpid()) + '.log'):
                        output = worker.Apply(data.chrom, varnames, sorted(varpos), data.getFamSamples(item))
                    if len(output) == 0:
                        # core algorithm failed
                        with env.lock:
                            env.chperror_counter.value += 1
                    else:
                        for line in output:
                            # line: [fid, sid, hap1, hap2]
                            data[line[1]] = (line[2], line[3])
                            if len(line[2]) > data.superMarkerCount:
                                data.superMarkerCount = len(line[2])
                        with env.lock:
                            env.mendelerror_counter.value += worker.countMendelianErrors()
                            env.recomb_counter.value += worker.countRecombs()
                    # FIXME: it is SWIG's (2.0.11) fault not to properly distroy the object "worker"
                    # So there is a memory leak here which I tried to partially handle on C++
            except Exception as e:
                return -1
        # Reformat sample genotypes 
        for person in data:
            if type(data[person]) is not tuple:
                data[person] = self.missings
            if data.superMarkerCount >= 2:
                diff = data.superMarkerCount - len(data[person][0])
                data[person] = [(x, y) for x, y in \
                                zip(data[person][0] + self.missing * diff, data[person][1] + self.missing * diff)]
        return 0


class LinkageWriter:
    def __init__(self):
        self.chrom = self.prev_chrom = self.name = self.distance = None
        self.reset()

    def apply(self, data):
        # input data receieved, call it a "success"
        with env.lock:
            env.success_counter.value += 1
        if self.chrom != self.prev_chrom:
            if self.prev_chrom is None:
                self.prev_chrom = self.chrom
            else:
                # new chrom entered,
                # commit whatever is in buffer before accepting new data
                self.commit()
        # write output
        position = str(data.getMidPosition())
        if data.superMarkerCount <= 1:
            gs = [data[s] for s in data.samples]
            if len(set(gs)) == 1:
                # everyone is missing (or monomorphic)
                return 2
            self.output += env.delimiter.join([self.chrom, self.name, self.distance, position]) + \
              env.delimiter + env.delimiter.join(chain(*gs)) + "\n" 
        else:
            # have to expand each region into mutiple sub-regions to account for different recomb points
            gs = zip(*[data[s] for s in data.samples])
            idx = 0
            for g in gs:
                if len(set(g)) == 1:
                    continue
                idx += 1
                self.output += \
                  env.delimiter.join([self.chrom, '{}[{}]'.format(self.name, idx), self.distance, position]) + \
                  env.delimiter + env.delimiter.join(chain(*g)) + "\n"
            if idx == 0:
                return 2
        if self.counter < env.batch:
            self.counter += data.superMarkerCount
        else:
            self.commit()
        return 0

    def commit(self):
        if not self.output:
            return
        with env.lock:
            with open(os.path.join(env.cache_dir, '{}.chr{}.tped'.format(env.output, self.prev_chrom)),
                      'a') as f:
                f.write(self.output)
        self.reset()

    def reset(self):
        self.output = ''
        self.counter = 0
        self.prev_chrom = self.chrom
            
    def getRegion(self, region):
        self.chrom = region[0]
        self.name, self.distance = region[3:]


class EncoderWorker(Process):
    def __init__(self, wid, queue, data, extractor, coder, writer):
        Process.__init__(self)
        self.worker_id = wid
        self.queue = queue
        self.data = data
        self.extractor = extractor
        self.coder = coder
        self.writer = writer

    def report(self):
        env.log('Processing {:,d} units with {:,d} super markers being generated ...'.\
                format(env.total_counter.value, env.success_counter.value), flush = True)
        
    def run(self):
        while True:
            try:
                region = self.queue.get()
                if region is None:
                    self.writer.commit()
                    self.report()
                    break
                else:
                    with env.lock:
                        env.total_counter.value += 1
                    self.extractor.getRegion(region)
                    self.writer.getRegion(region)
                    for m in [self.extractor, self.coder, self.writer]:
                        status = m.apply(self.data)
                        with env.lock:
                            if status == -1:
                                # previous module failed
                                env.chperror_counter.value += 1
                            if status == 1:
                                env.null_counter.value += 1
                            if status == 2:
                                env.trivial_counter.value += 1
                        if status != 0:
                            break
                    if env.total_counter.value % (env.batch * env.jobs) == 0:
                        self.report()
            except KeyboardInterrupt:
                break


def main(args):
    '''the main encoder function'''
    checkParams(args)
    # FIXME: add Di's resources & MAC version resource
    downloadResources([('http://tigerwang.org/uploads/genemap.txt', env.resource_dir),
                       ('http://tigerwang.org/uploads/transpose.pl', env.resource_bin),
                       ('http://tigerwang.org/uploads/runMlink.pl', env.resource_bin)])
    cache = Cache(env.cache_dir, env.output)
    # STEP 1: write encoded data to TPED format
    if not args.vanilla and cache.check():
        env.log('Loading data from archive ...')
        cache.load()
    else:
        # load VCF file header
        checkVCFBundle(args.vcf)
        cache.clear(pre = env.output, ext =".tped")
        try:
            vs = cstatgen.VCFstream(args.vcf)
        except Exception as e:
            env.error("{}".format(e), exit = True)
        samples_vcf = vs.GetSampleNames()
        if len(samples_vcf) == 0:
            env.error("Fail to extract samples from [{}]".format(vcf), exit = True)
        env.log('{:,d} samples found in [{}]'.format(len(samples_vcf), args.vcf))
        checkSamples(samples_vcf, getColumn(args.tfam, 2))
        # load sample info from tfam file
        tfam = TFAMParser(args.tfam)
        if len(tfam.families) == 0:
            env.error('No family found in [{}]!'.format(args.tfam), exit = True)
        else:
            # samples and families have to be in both tfam file and vcf file
            samples = {k : tfam.samples[k] for k in tfam.samples if k in samples_vcf}
            families = {k : [x for x in tfam.families[k] if x in samples] for k in tfam.families}
            for k in families.keys():
                if len(families[k]) == 0:
                    del families[k]
            #
            if len(samples) == 0:
                env.error('No valid sample to process. '\
                          'Samples have to be in families, and present in both TFAM and VCF files.', exit = True)
            else:
                env.log('{:,d} families with a total of {:,d} samples will be processed'.\
                    format(len(tfam.families), len(samples)))
        rewriteFamfile(env.outputfam, samples, [x for x in samples_vcf if x in samples])
        # load blueprint
        with open(args.blueprint, 'r') as f:
            regions = [x.strip().split() for x in f.readlines()]
        env.log('Scanning for {:,d} pre-defined units in VCF file:'.format(len(regions)))
        env.jobs = max(min(args.jobs, len(regions)), 1)
        regions.extend([None] * env.jobs)
        queue = Queue()
        try:
            faulthandler.enable(file=open(env.tmp_log + '.SEGV', 'w'))
            for i in regions:
                queue.put(i)
            jobs = [EncoderWorker(
                i,
                queue,
                RData(samples, families, samples_vcf),
                RegionExtractor(args.vcf),
                MarkerMaker(args.size),
                LinkageWriter()
                ) for i in range(env.jobs)]
            for j in jobs:
                j.start()
            for j in jobs:
                j.join()
            faulthandler.disable()
        except KeyboardInterrupt:
            # FIXME: need to properly close all jobs
            raise ValueError("Use 'killall {}' to properly terminate all processes!".format(env.prog))
        else:
            env.log('{:,d} variants processed with {:,d} super markers generated; '\
                '{:,d} Mendelian inconsistencies and {:,d} recombination events handled\n'.\
                format(env.variants_counter.value,
                       env.success_counter.value,
                       env.mendelerror_counter.value,
                       env.recomb_counter.value), flush = True)
            if env.triallelic_counter.value:
                env.log('{:,d} tri-allelic loci were ignored'.format(env.triallelic_counter.value))
            if env.null_counter.value:
                env.log('{:,d} units ignored due to absence in VCF file'.format(env.null_counter.value))
            if env.trivial_counter.value:
                env.log('{:,d} units ignored due to absence of variants in samples'.format(env.trivial_counter.value))
            fatal_errors = 0
            try:
                # Error msg from C++ extension
                os.system("cat {}/*.* > {}".format(env.tmp_dir, env.tmp_log))
                fatal_errors = wordCount(env.tmp_log)['fatal']
            except KeyError:
                pass
            if env.chperror_counter.value:
                env.error("{:,d} super markers failed to be generated due to genotype coding problems!".\
                          format(env.chperror_counter.value))
            if fatal_errors:
                env.error("{:,d} or more super markers failed to be generated due to runtime errors!".\
                          format(fatal_errors))
            env.log('Archiving to directory [{}]'.format(env.cache_dir))
            cache.write(pre = env.output, ext = '.tped', files = [env.outputfam],
                        otherfiles = [args.vcf, args.tfam, args.blueprint])
    env.jobs = args.jobs
    # STEP 2: write to PLINK or mega2 format
    tpeds = [os.path.join(env.cache_dir, item) for item in os.listdir(env.cache_dir) if item.startswith(env.output) and item.endswith('.tped')]
    if 'plink' in args.format:
        env.log('Saving data to directory [PLINK] ...')
        formatPlink(tpeds, [env.outputfam] * len(tpeds), 'PLINK')
    if 'mega2' in args.format:
        env.log('Saving data to directory [MEGA2] ...')
        formatMega2('MEGA2/{}*'.format(env.output))
    if 'mlink' in args.format:
        env.log('Saving data to directory [MLINK] ...')
        formatMlink(tpeds, [env.outputfam] * len(tpeds), 'MLINK')
    runMlink()
    plotMlink()
    #please implement this section
    # STEP final: clean up unwanted files
    #for item in args.format:
    #    removeFiles(item.upper(), exclude = env.formats[item])
    #removeFiles(env.cache_dir, exclude = ['.cache'])
