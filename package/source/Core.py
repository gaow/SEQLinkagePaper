#!/usr/bin/python
# Copyright (c) 2013, Gao Wang <gaow@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

from SEQLinco.Utils import *
from multiprocessing import Process, Queue, Lock, Value
from collections import Counter
import sys, gc

if sys.version_info.major == 2:
    from SEQLinco.libmped import chp_py2 as MM
else:
    from SEQLinco.libmped import chp_py3 as MM

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

def rewriteFamfile(tfam, samples):
    os.rename(tfam, tfam + '.bak')
    with open(tfam, 'w') as f:
        f.write('\n'.join(sorted(['\t'.join(x) for x in samples.values()])))
        f.write('\n')

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
        fams = {}
        observedFounders = {}
        expectedFounders = {}
        samples = {}
        with open(tfam, 'r') as f:
            for line in f.readlines():
                line = line.split()
                # collect sample info
                samples[line[1]] = [line[0], line[1], line[2], line[3], line[4], line[5]]
                # collect family member
                self.__add_or_app(fams, line[0], line[1])
                # collect founders for family
                if line[2] in env.ped_missing or line[3] in env.ped_missing:
                    self.__add_or_app(observedFounders, line[0], line[1])
                else:
                    self.__add_or_app(expectedFounders, line[0], (line[1], line[2], line[3]))
        # if a sample's expected founder is not observed in tfam
        # then the sample itself is a founder 
        for k in expectedFounders.keys():
            if k not in observedFounders:
                for item in expectedFounders[k]:
                    samples[item[0]][2] = samples[item[0]][3] = "0"
                observedFounders[k] = [item[0] for item in expectedFounders[k]]
                continue
            for item in expectedFounders[k]:
                if item[1] not in observedFounders[k] \
                  and item[2] not in observedFounders[k]:
                    samples[item[0]][2] = samples[item[0]][3] = "0"
                    observedFounders[k].append(item[0])
        # now remove trivial families 
        # 1. the family does not have observedFounders
        # 2. families having founders only, no offspring, i.e., observedFounder is all samples in the fam
        for k in fams.keys():
            if k not in observedFounders:
                del fams[k]
                continue
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
    def __init__(self, samples, families, sample_names):
        # a dict of {sid:[fid, pid, mid, sex, trait], ...}
        self.samples = samples
        # a list of sample names matching the order of apparence in VCF file
        self.sample_names = sample_names
        # a dict of {fid:[member names], ...}
        # reorder family samples based on order of VCF file
        self.families = {}
        self.famsampidx = {}
        self.famvaridx = {}
        for k in families:
            self.families[k] = [x for x in self.sample_names if x in families[k]]
            self.famsampidx[k] = [i for i, x in enumerate(self.sample_names) if x in families[k]]
        self.reset()

    def reset(self):
        for item in self.sample_names:
            self[item] = []
        self.variants = []
        for k in self.families:
            self.famvaridx[k] = []

    def getMidPosition(self):
        if len(self.variants) == 0:
            return None
        return sum([int(x[1]) for x in self.variants]) / len(self.variants)

    def getFamVariants(self, fam, style = None):
        famvar = [item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam]]
        if style is None:
            return famvar
        elif style == "map":
            names = []
            pos = []
            for item in famvar:
                names.append("{}-{}".format(item[-1], item[1]))
                pos.append(item[1])
            return names, pos
        else:
            raise ValueError("Unknown style '{}'".format(style))

    def getFamSamples(self, fam):
        output = [[]] * len(self.families[fam])
        for idx, item in enumerate(self.families[fam]):
            output[idx] = self.samples[item][:-1]
            output[idx].extend(self[item])
        return output


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
        self.lock = Lock()
        
    def apply(self, data):
        # Clean up
        data.reset()
        lines = [self.parseVCFline(x, exclude = self.exclude_idx) for x in \
                  runCommand('tabix {} {}:{}-{}'.\
                             format(self.vcf, self.chrom, self.startpos, self.endpos)).strip().split('\n')]
        lines = filter(None, lines)
        if len(lines) == 0:
            return False
        # check if the first line's sample name matches with expected sample name
        if not len(lines[0][1]) == len(data.sample_names):
            raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name, len(lines[0][1]), len(data.sample_names)))
        #
        # FIXME: codes below do not read efficient
        #
        # for each variant
        for var_id, line in enumerate(lines):
            data.variants.append(line[0])
            # for each family assign member genotype if the site is non-trivial in the fam 
            for k in data.families:
                gs = [y for idx, y in enumerate(line[1]) if idx in data.famsampidx[k]]
                if len(set(''.join([x for x in gs if x != "00"]))) == 1:
                    # skip monomorphic gs
                    continue
                else:
                    data.famvaridx[k].append(var_id)
                    for person, g in zip(data.families[k], gs):
                        data[person].append(g)
        return True

            
    def getRegion(self, region):
        self.chrom, self.startpos, self.endpos, self.name, self.distance = region
        if self.chrom.lower() in ['x','23','chrx','chr23']:
            if self.xchecker.check(int(self.startpos)) or self.xchecker.check(int(self.endpos)): 
                self.chrom = 'XY'
        if self.chrom.lower() in ['y','24','chry','chr24']:
            if self.ychecker.check(int(self.startpos)) or self.ychecker.check(int(self.endpos)):
                self.chrom = 'XY'

    
    def parseVCFline(self, line, exclude = []):
        if len(line) == 0:
            return None
        line = line.split('\t')
        # Skip tri-allelic variant
        if "," in line[4]:
            with self.lock:
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
        return (line[0], line[1], line[3], line[4], self.name), gs

    
class MarkerMaker:
    def __init__(self, wsize):
        self.size = wsize
        self.lock = Lock()
        env.tmp_log = os.path.join(env.tmp_dir, env.output)
        # adjust physical distance to map distance, 1 / 1 million * 10 (times 10 to fake them apart)
        self.position_adj = 1 / float(100000)
        self.missing = env.delimiter.join(["0", "0"])

    def apply(self, data):
        for item in data.families:
            try:
                vnames, vpos = data.getFamVariants(item, style = "map")
                if len(vnames) < 2:
                    # no more than 2 variants found in family
                    for person in data.families[item]:
                        if len(data[person]) == 0:
                            data[person] = self.missing
                        else:
                            data[person] = env.delimiter.join(list(data[person][0]))
                else:
                    worker = MM.CHP(self.size, self.position_adj)
                    with stdoutRedirect(to = env.tmp_log + str(os.getpid()) + '.log'):
                        output = worker.Apply(data.variants[0][0],
                                              vnames,
                                              sorted(vpos),
                                              data.getFamSamples(item))
                    if len(output) == 0:
                        with self.lock:
                            env.chperror_counter.value += 1
                    else:
                        for line in output:
                            # line: [fid, sid, hap1, hap2]
                            data[line[1]] = env.delimiter.join([line[2], line[3]])
                        with self.lock:
                            env.mendelerror_counter.value += worker.countMendelianErrors()
                            env.recomb_counter.value += worker.countRecombs()
            except Exception as e:
                raise
                return False
        # Fill in skipped samples due to some errors 
        for person in data:
            if type(data[person]) is not str:
                data[person] = self.missing
        return True


class LinkageWriter:
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
        # FIXME: currently only applies to one marker per region theme
        self.output += env.delimiter.join(
            [self.chrom, self.name, self.distance, str(data.getMidPosition())] + \
                                          [data[s] for s in data.sample_names]) + '\n'
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
    def __init__(self, wid, queue, data, extractor, coder, writer):
        Process.__init__(self)
        self.worker_id = wid
        self.queue = queue
        self.data = data
        self.extractor = extractor
        self.coder = coder
        self.writer = writer
        self.lock = Lock()

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
                    with self.lock:
                        env.total_counter.value += 1
                    self.extractor.getRegion(region)
                    self.writer.getRegion(region)
                    for m in [self.extractor, self.coder, self.writer]:
                        if not m.apply(self.data):
                            # previous module failed
                            break
                    if env.total_counter.value % (env.batch * env.jobs) == 0:
                        self.report()
            except KeyboardInterrupt:
                break


def main(args):
    '''the main encoder function'''
    status = checkParams(args)
    # FIXME: add Di's resources & MAC version resource
    status = downloadResources([('http://tigerwang.org/uploads/genemap.txt', env.resource_dir),
                               ('http://tigerwang.org/uploads/tabix', env.resource_bin),
                               ('http://tigerwang.org/uploads/bgzip', env.resource_bin)])
    cache = Cache(env.cache_dir, env.output)
    # STEP 1: write encoded data to TPED format
    if not args.vanilla and cache.check():
        env.log('Loading data from archive ...')
        cache.load()
    else:
        args.vcf = indexVCF(args.vcf)
        samples_vcf = extractSamplenames(args.vcf)
        status = checkSamples(samples_vcf, getColumn(args.tfam, 2))
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
        #
        status = rewriteFamfile(args.tfam, samples)
        with open(args.blueprint, 'r') as f:
            regions = [x.strip().split() for x in f.readlines()]
        env.log('Scanning [{}] for {:,d} pre-defined units'.format(args.vcf, len(regions)))
        env.jobs = max(min(args.jobs, len(regions)), 1)
        regions.extend([None] * env.jobs)
        queue = Queue()
        try:
            for i in regions:
                queue.put(i)
            jobs = [EncoderWorker(
                i,
                queue,
                RData(samples, families, sample_names = [x for x in samples_vcf if x in samples]),
                RegionExtractor(args.vcf, [idx + 9 for idx, x in enumerate(samples_vcf) if x not in samples]),
                MarkerMaker(args.size),
                LinkageWriter()
                ) for i in range(env.jobs)]
            for j in jobs:
                j.start()
            for j in jobs:
                j.join()
        except KeyboardInterrupt:
            # FIXME: need to properly close all jobs
            raise ValueError("Use 'killall {}' to properly terminate all processes!".format(env.prog))
        else:
            env.log('{:,d} units processed with {:,d} super markers generated; '\
                '{:,d} tri-allelic loci ignored, {:,d} Mendelian errors fixed, '\
                '{:,d} recombination events detected.\n'.\
                format(env.total_counter.value,
                       env.success_counter.value,
                       env.triallelic_counter.value,
                       env.mendelerror_counter.value,
                       env.recomb_counter.value), flush = True)
            if env.chperror_counter.value > 0:
                env.error("{:,d} runtime error captured!".format(env.chperror_counter.value))
            env.log('Archiving to directory [{}]'.format(env.cache_dir))
            cache.write(pre = env.output, ext = '.tped', otherfiles = [args.vcf, args.tfam, args.blueprint])
    # STEP 2: write to PLINK or mega2 format
    sys.exit()
    env.jobs = args.jobs
    tpeds = [os.path.join(env.cache_dir, item) for item in os.listdir(env.cache_dir) if item.startswith(env.output) and item.endswith('.tped')]
    if 'plink' in args.format:
        env.log('Saving data to directory [PLINK] ...')
        formatPlink(tpeds, [args.tfam] * len(tpeds), 'PLINK')
    if 'mega2' in args.format:
        env.log('Saving data to directory [MEGA2] ...')
        formatMega2('MEGA2/{}*'.format(env.output))
    if 'mlink' in args.format:
        env.log('Saving data to directory [MLINK] ...')
        formatMlink(tpeds, [args.tfam] * len(tpeds), 'MLINK')
    runMlink()
    plotMlink()
    #please implement this section
    # STEP final: clean up unwanted files
    #for item in args.format:
    #    removeFiles(item.upper(), exclude = env.formats[item])
    #removeFiles(env.cache_dir, exclude = ['.cache'])
