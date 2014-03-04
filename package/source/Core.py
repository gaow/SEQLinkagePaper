#!/usr/bin/python
# Copyright (c) 2013, Gao Wang <gaow@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
from __future__ import print_function
from SEQLinco.Utils import *
from SEQLinco.Runner import *
from multiprocessing import Process, Queue
from collections import Counter, OrderedDict, defaultdict
import itertools
from copy import deepcopy
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
        env.log('{:,d} samples found in TFAM file but not in VCF file:\n{}'.\
                           format(len(b_not_a), '\n'.join(b_not_a)))
    if a_not_b:
        env.log('{:,d} samples in VCF file will be ignored due to absence in TFAM file'.format(len(a_not_b)))
    return a_not_b, b_not_a

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

    def write(self, pres = [], exts = [], files = [], otherfiles = []):
        '''Add files to cache'''
        if not self.check():
            with tarfile.open(self.cache_name, 'w:bz2') as f:
                for item in os.listdir(self.cache_dir):
                    if ((any([item.endswith(x) for x in exts]) or len(exts) == 0) \
                        and (any([item.startswith(x) for x in pres]) or len(pres) == 0)) \
                      or item in files:
                        f.add(os.path.join(self.cache_dir, item), arcname=item)
            otherfiles.append(self.cache_name)
            signatures = ['{}\t{}'.format(x, calculateFileMD5(x)) for x in otherfiles]
            with open(self.cache_info, 'w') as f:
                f.write('\n'.join(signatures))

    def clear(self, pres = [], exts = []):
        for fl in [self.cache_info, self.cache_name]:
            if os.path.isfile(fl):
                os.remove(fl)
        #
        for pre, ext in itertools.product(pres, exts): 
            for fl in glob.glob(os.path.join(self.cache_dir, pre) +  "*" + ext): 
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
    def __init__(self, tfam = None):
        self.families, self.samples, self.graph = self.__parse(tfam)
        self.families_sorted = OrderedDict([(k,[]) for k in self.families])

    def is_founder(self, sid):
        return self.samples[sid][2] == "0" and self.samples[sid][3] == "0"

    def get_parents(self, sid):
        return self.samples[sid][2], self.samples[sid][3]

    def add_member(self, info):
        '''member is one line of TFAM file, [fid, sid, pid, mid, sex, pheno]'''
        self.samples[info[1]] = info
        self.families_sorted[info[0]] = []
        if info[0] not in self.families:
            self.families[info[0]] = [info[1]] 
        else:
            if info[1] not in self.families[info[0]]:
                self.families[info[0]].append(info[1])
        self.__update_graph(self.graph, info)

    def get_member_idx(self, sid):
        '''integer index, reflecting the order of the sample collected'''
        return self.samples.keys().index(sid)

    def get_members(self):
        return self.samples.keys()

    def sort_family(self, famid):
        '''sort samples in family such that founders precede non-founders'''
        if not self.families_sorted[famid]:
            self.families_sorted[famid] = self.__kahn_sort(famid)
        assert sorted(self.families_sorted[famid]) == sorted(self.families[famid])
        return self.families_sorted[famid]
            
    def __kahn_sort(self, famid):
        '''algorithm first described by Kahn (1962); implemented by Di Zhang'''
        sorted_names = []
        S_no_parents = filter(lambda x: True if self.is_founder(x) else False, self.families[famid])
        graph = self.graph[famid].copy()
        while(S_no_parents):
            n = S_no_parents.pop()
            sorted_names.append(n)
            if n not in graph:
                continue
            offsprings = graph.pop(n)
            for m in offsprings:
                father, mother = self.get_parents(m)
                if father not in graph and mother not in graph:
                    S_no_parents.append(m)
        if graph:
            raise ValueError("There is a loop in the pedigree: {}\n".format(' '.join(graph.keys())))
        else:
            return sorted_names

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

    def __update_graph(self, g, info):
        if info[2] != "0" and info[3] != "0":
            g[info[0]][info[2]].append(info[1]) 
            g[info[0]][info[3]].append(info[1]) 

    def __parse(self, tfam):
        '''Rules:
        1. samples have to have unique names
        2. both parents for a non-founder should be available
        3. founders should have at least one offspring'''
        fams = OrderedDict()
        samples = OrderedDict()
        graph = defaultdict(lambda : defaultdict(list))
        if tfam is None:
            return fams, samples, graph
        observedFounders = {}
        expectedParents = {}
        #
        # Load TFAM file
        #
        with open(tfam, 'r') as f:
            for idx, line in enumerate(f.readlines()):
                line = line.split()
                if len(line) != 6:
                    env.error("skipped line {} (has {} != 6 columns!)".format(idx, len(line)))
                    continue
                if line[1] in samples:
                    env.error("skipped line {} (duplicate sample name '{}' found!)".format(idx, line[1]))
                    continue
                # collect sample line 
                samples[line[1]] = [line[0], line[1], line[2], line[3], line[4], line[5]]
                # collect family member
                self.__add_or_app(fams, line[0], line[1])
                # collect founders for family
                if line[2] in env.ped_missing and line[3] in env.ped_missing:
                    self.__add_or_app(observedFounders, line[0], line[1])
                else:
                    self.__add_or_app(expectedParents, line[0], (line[1], line[2], line[3]))
        #
        # Check sample parents
        #
        for k in expectedParents:
            for person in expectedParents[k]:
                if not (person[1] in fams[k] and person[2] in fams[k]):
                    env.error("Cannot find parents ({} and {}) of {} in [{}]!".\
                              format(person[1], person[2], person[0], tfam), exit = True)
                    # missing both parents, make it a founder
                    # samples[person[0]][2] = samples[person[0]][3] = "0"
                    # observedFounders[k].append(person[0])
                if person[1] in fams[k] and not person[2] in fams[k]:
                    env.error("Cannot find mother ({}) of {} in [{}]!".\
                              format(person[2], person[0], tfam), exit = True)
                    # missing mother, mask as zero
                    # samples[person[0]][3] = "0"
                if not person[1] in fams[k] and person[2] in fams[k]:
                    env.error("Cannot find father ({}) of {} in [{}]!".\
                              format(person[1], person[0], tfam), exit = True)
                    # missing father, mask as zero
                    # samples[person[0]][2] = "0"
        #
        # Remove trivial families 
        #
        for k in fams.keys():
            if Counter(observedFounders[k]) == Counter(fams[k]):
                del fams[k]
                continue
        #
        valid_samples = []
        for value in fams.values():
            valid_samples.extend(value)
        samples = {k : samples[k] for k in valid_samples}
        #
        for item in samples.values():
            self.__update_graph(graph, item)
        return fams, samples, graph
    
    
class RData(dict):
    def __init__(self, samples_vcf, tfam):
        # tfam.samples: a dict of {sid:[fid, pid, mid, sex, trait], ...}
        # tfam.families: a dict of {fid:[s1, s2 ...], ...}
        self.tfam = tfam
        # samples have to be in both vcf and tfam data
        self.samples = OrderedDict([(k, tfam.samples[k]) for k in samples_vcf if k in tfam.samples])
        # a dict of {fid:[member names], ...}
        self.families = {k : [x for x in self.samples if x in tfam.families[k]] for k in tfam.families}
        # a dict of {fid:[idx ...], ...}
        self.famsampidx = {}
        # a dict of {fid:[maf1, maf2 ...]}
        self.maf = OrderedDict()
        # reorder family samples based on order of VCF file
        for k in self.families.keys():
            if len(self.families[k]) == 0:
                # skip families having no samples in VCF file
                del self.families[k]
            else:
                self.famsampidx[k] = [i for i, x in enumerate(samples_vcf) if x in self.families[k]]
        # a dict of {fid:[idx ...], ...}
        self.famvaridx = {}
        self.reset()

    def reset(self):
        for item in self.samples:
            self[item] = []
        self.variants = []
        self.chrom = None
        for k in self.families:
            self.famvaridx[k] = []
        self.maf = OrderedDict() 
        # superMarkerCount is the max num. of recombinant fragments among all fams   
        self.superMarkerCount = 0

    def getMidPosition(self):
        if len(self.variants) == 0:
            return None
        return sum([x[1] for x in self.variants]) / len(self.variants)

    def getFamVariants(self, fam, style = None):
        if style is None:
            return [item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam]]
        elif style == "map":
            names = []
            pos = []
            mafs = []
            for idx in self.famvaridx[fam]:
                names.append("V{}".format(idx))
                pos.append(self.variants[idx][1])
                mafs.append(self.variants[idx][-1])
            return names, pos, mafs
        else:
            raise ValueError("Unknown style '{}'".format(style))

    def getFamSamples(self, fam):
        nvar = len([item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam]])
        output = [[]] * len(self.tfam.families[fam])
        for idx, item in enumerate(self.tfam.sort_family(fam)):
            # sample info, first 5 columns of ped
            output[idx] = self.tfam.samples[item][:-1]
            # sample genotypes
            if item in self.samples:
                output[idx].extend(self[item])
            else:
                output[idx].extend(["00"] * nvar)
        return output

    
class RegionExtractor:
    '''Extract given genomic region from VCF
    converting genotypes into dictionary of
    genotype list'''
    def __init__(self, filename, build = env.build, chr_prefix = None, allele_freq_info = None):
        self.vcf = cstatgen.VCFstream(filename)
        self.chrom = self.startpos = self.endpos = self.name = self.distance = None
        self.chr_prefix = chr_prefix
        # name of allele frequency meta info
        self.af_info = allele_freq_info
        self.xchecker = PseudoAutoRegion('X', build)
        self.ychecker = PseudoAutoRegion('Y', build)
        
    def apply(self, data):
        # Clean up
        data.reset()
        data.chrom = self.chrom
        self.vcf.Extract(self.chrom, self.startpos, self.endpos)
        varIdx = 0
        # for each variant site
        while (self.vcf.Next()):
            # skip tri-allelic sites
            if not self.vcf.IsBiAllelic():
                with env.triallelic_counter.get_lock():
                    env.triallelic_counter.value += 1
                continue
            # check if the line's sample number matches the entire VCF sample number 
            if not self.vcf.CountSampleGenotypes() == self.vcf.sampleCount:
                raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name, self.vcf.CountSampleGenotypes(), self.vcf.sampleCount))
            # valid line found, get variant info
            try:
                maf = float(self.vcf.GetInfo(self.af_info)) if self.af_info else None 
                if maf > 0.5:
                    maf = 1 - maf
            except Exception as e:
                raise ValueError("VCF line {}:{} does not have valid allele frequency field {}!".\
                                 format(self.vcf.GetChrom(), self.vcf.GetPosition(), self.af_info))
            data.variants.append([self.vcf.GetChrom(), self.vcf.GetPosition(), self.name, maf])
            # for each family assign member genotype if the site is non-trivial to the family
            for k in data.families:
                gs = self.vcf.GetGenotypes(data.famsampidx[k])
                if len(set(''.join([x for x in gs if x != "00"]))) <= 1:
                    # skip monomorphic gs
                    continue
                else:
                    # this variant is found in the family
                    data.famvaridx[k].append(varIdx)
                    for person, g in zip(data.families[k], gs):
                        data[person].append(g)
            varIdx += 1
        # 
        if varIdx == 0:
            return 1
        else:
            with env.variants_counter.get_lock():
                env.variants_counter.value += varIdx 
            return 0
            

    def getRegion(self, region):
        self.chrom, self.startpos, self.endpos, self.name, self.distance = region
        self.startpos = int(self.startpos)
        self.endpos = int(self.endpos)
        if self.chrom in ['X','23']:
            if self.xchecker.check(self.startpos) or self.xchecker.check(self.endpos): 
                self.chrom = 'XY'
        if self.chrom in ['Y','24']:
            if self.ychecker.check(self.startpos) or self.ychecker.check(self.endpos):
                self.chrom = 'XY'
        if self.chr_prefix and not self.chrom.startswith(self.chr_prefix):
            self.chrom = self.chr_prefix + self.chrom

    
class MarkerMaker:
    def __init__(self, wsize):
        self.missings = ("0", "0")
        self.missing = "0"
        self.gtconv = {'1':0, '2':1}
        self.haplotyper = cstatgen.HaplotypingEngine(verbose = env.debug)
        self.coder = cstatgen.HaplotypeCoder(wsize)

    def apply(self, data):
        # temp raw haplotype, maf and variant names data
        haplotypes = OrderedDict() 
        mafs = {}
        varnames = {}
        try:
            # haplotyping plus collect found allele counts
            # and computer founder MAFS
            self.__Haplotype(data, haplotypes, mafs, varnames)
            # recoding the genotype of the region
            self.__CodeHaplotypes(data, haplotypes, mafs, varnames)
        except Exception as e:
            raise
            return -1
        self.__FormatHaplotypes(data)
        return 0
    
    def __Haplotype(self, data, haplotypes, mafs, varnames):
        '''genetic haplotyping. haplotypes stores per family data'''
        # FIXME: it is SWIG's (2.0.12) fault not to properly destroy the object "Pedigree" in "Execute()"
        # So there is a memory leak here which I tried to partially handle on C++
        #
        # Per family haplotyping
        #
        for item in data.families:
            varnames[item], positions, vcf_mafs = data.getFamVariants(item, style = "map")
            if len(varnames[item]) == 0:
                for person in data.families[item]:
                    data[person] = self.missings
                continue
            if env.debug:
                with env.lock:
                    sys.stderr.write('\n'.join(['\t'.join(x) for x in data.getFamSamples(item)]) + '\n\n')
            # haplotyping
            with env.lock:
                with stdoutRedirect(to = env.tmp_log + str(os.getpid()) + '.log'):
                    haplotypes[item] = self.haplotyper.Execute(data.chrom, varnames[item],
                                                               sorted(positions), data.getFamSamples(item))[0]
            if len(haplotypes[item]) == 0:
                # C++ haplotyping implementation failed
                with env.chperror_counter.get_lock():
                    env.chperror_counter.value += 1
            # either use privided MAF or computer MAF
            if all(vcf_mafs):
                for idx, v in enumerate(varnames[item]):
                    if v not in mafs:
                        mafs[v] = vcf_mafs[idx]
            else:
                # count founder alleles
                for hap in haplotypes[item]:
                    if not data.tfam.is_founder(hap[1]):
                        continue
                    for idxv, v in enumerate(varnames[item]):
                        if v not in mafs:
                            # [#alt, #haplotypes]
                            mafs[v] = [0, 0]
                        gt = hap[2 + idxv][1] if hap[2 + idxv][0].isupper() else hap[2 + idxv][0]
                        if not gt == "?":
                            mafs[v][0] += self.gtconv[gt]
                            mafs[v][1] += 1.0 
        #
        # Compute founder MAFs
        #
        for v in mafs:
            if type(mafs[v]) is not list:
                continue
            mafs[v] = mafs[v][0] / mafs[v][1] if mafs[v][1] > 0 else 0.0
        if env.debug:
            with env.lock:
                print("variant mafs = ", mafs, "\n", file = sys.stderr)

    def __CodeHaplotypes(self, data, haplotypes, mafs, varnames):
        self.coder.Execute(haplotypes.values(), [[mafs[v] for v in varnames[item]] for item in haplotypes])
        if env.debug:
            with env.lock:
                self.coder.Print()
        # line: [fid, sid, hap1, hap2]
        for line in self.coder.GetHaplotypes():
            if not line[1] in data:
                # this sample is not in VCF file. Every variant site should be missing
                # they have to be skipped for now
                continue 
            data[line[1]] = (line[2], line[3])
            if len(line[2]) > data.superMarkerCount:
                data.superMarkerCount = len(line[2])
        # get MAF
        for item in haplotypes:
            data.maf[item] = self.coder.GetAlleleFrequencies(item)
        if env.debug:
            with env.lock:
                print("marker freqs = ", data.maf, "\n", file = sys.stderr)
    
    def __FormatHaplotypes(self, data):
        # Reformat sample genotypes 
        for person in data:
            if type(data[person]) is not tuple:
                data[person] = self.missings
            if data.superMarkerCount >= 2:
                diff = data.superMarkerCount - len(data[person][0])
                data[person] = [(x, y) for x, y in \
                                zip(data[person][0] + self.missing * diff, data[person][1] + self.missing * diff)]
        

class LinkageWriter:
    def __init__(self, num_missing_append = 0):
        self.chrom = self.prev_chrom = self.name = self.distance = None
        self.reset()
        self.missings = ["0", "0"]
        self.num_missing = num_missing_append

    def apply(self, data):
        if self.chrom != self.prev_chrom:
            if self.prev_chrom is None:
                self.prev_chrom = self.chrom
            else:
                # new chrom entered,
                # commit whatever is in buffer before accepting new data
                self.commit()
        # write tped output
        position = str(data.getMidPosition())
        if data.superMarkerCount <= 1:
            gs = [data[s] for s in data.samples]
            if len(set(gs)) == 1:
                # everyone's genotype is the same (most likely missing or monomorphic)
                return 2
            self.tped += env.delimiter.join([self.chrom, self.name, self.distance, position] + \
                list(itertools.chain(*gs)) + self.missings*self.num_missing) + "\n"
        else:
            # have to expand each region into mutiple sub-regions to account for different recomb points
            gs = zip(*[data[s] for s in data.samples])
            idx = 0
            for g in gs:
                if len(set(g)) == 1:
                    continue
                idx += 1
                self.tped += \
                  env.delimiter.join([self.chrom, '{}[{}]'.format(self.name, idx), self.distance, position] + \
                  list(itertools.chain(*g)) + self.missings*self.num_missing) + "\n"
            if idx == 0:
                # everyone's genotype is the same (most likely missing or monomorphic)
                return 2
        # write freq file output
        if data.superMarkerCount > 1:
            for k in data.maf:
                for idx in range(data.superMarkerCount):
                    if idx < len(data.maf[k]):
                        self.freq += env.delimiter.join([k, '{}[{}]'.format(self.name, idx + 1)] + \
                                                        map(str, data.maf[k][idx])) + "\n" 
        else:
            for k in data.maf:
                self.freq += env.delimiter.join([k, self.name] + map(str, data.maf[k][0])) + "\n" 
        if self.counter < env.batch:
            self.counter += data.superMarkerCount
        else:
            self.commit()
        return 0

    def commit(self):
        if self.tped:
            with env.lock:
                with open(os.path.join(env.cache_dir, '{}.chr{}.tped'.format(env.output, self.prev_chrom)),
                          'a') as f:
                    f.write(self.tped)
        if self.freq:
            with env.lock:
                with open(os.path.join(env.cache_dir, '{}.chr{}.freq'.format(env.output, self.prev_chrom)),
                          'a') as f:
                    f.write(self.freq)
        self.reset()

    def reset(self):
        self.tped = ''
        self.freq = ''
        self.counter = 0
        self.prev_chrom = self.chrom
            
    def getRegion(self, region):
        self.chrom = region[0]
        self.name, self.distance = region[3:]


class EncoderWorker(Process):
    def __init__(self, queue, length, data, extractor, coder, writer):
        Process.__init__(self)
        self.queue = queue
        self.numGrps = float(length)
        self.data = data
        self.extractor = extractor
        self.maker = coder
        self.writer = writer

    def report(self):
        env.log('{:,d} units processed {{{:.2%}}} ...'.\
                format(env.success_counter.value, env.total_counter.value / self.numGrps), flush = True)
        
    def run(self):
        while True:
            try:
                region = self.queue.get()
                if region is None:
                    self.writer.commit()
                    self.report()
                    # total mendelian errors found
                    with env.mendelerror_counter.get_lock():
                        env.mendelerror_counter.value += self.maker.haplotyper.CountMendelianErrors()
                    # total recombination events found
                    with env.recomb_counter.get_lock():
                        env.recomb_counter.value += self.maker.coder.CountRecombinations()
                    break
                else:
                    with env.total_counter.get_lock():
                        env.total_counter.value += 1
                    self.extractor.getRegion(region)
                    self.writer.getRegion(region)
                    isSuccess = True
                    for m in [self.extractor, self.maker, self.writer]:
                        status = m.apply(self.data)
                        if status == -1:
                            with env.chperror_counter.get_lock():
                                # previous module failed
                                env.chperror_counter.value += 1
                        if status == 1:
                            with env.null_counter.get_lock():
                                env.null_counter.value += 1
                        if status == 2:
                            with env.trivial_counter.get_lock():
                                env.trivial_counter.value += 1
                        if status != 0:
                            isSuccess = False
                            break
                    if isSuccess:
                        with env.success_counter.get_lock():
                            env.success_counter.value += 1
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
        cache.clear(pres = [env.output], exts = [".tped", ".freq"])
        try:
            vs = cstatgen.VCFstream(args.vcf)
        except Exception as e:
            env.error("{}".format(e), exit = True)
        samples_vcf = vs.GetSampleNames()
        if len(samples_vcf) == 0:
            env.error("Fail to extract samples from [{}]".format(vcf), exit = True)
        env.log('{:,d} samples found in [{}]'.format(len(samples_vcf), args.vcf))
        samples_not_vcf = checkSamples(samples_vcf, getColumn(args.tfam, 2))[1]
        # load sample info 
        data = RData(samples_vcf, TFAMParser(args.tfam))
        if len(data.families) == 0:
            env.error('No valid family to process. ' \
                      'Families have to be at least trio with at least one member in VCF file.', exit = True)
        if len(data.samples) == 0:
            env.error('No valid sample to process. ' \
                      'Samples have to be in families, and present in both TFAM and VCF files.', exit = True)
        rewriteFamfile(env.outputfam, data.tfam.samples, data.samples.keys() + samples_not_vcf)
        # load blueprint
        with open(args.blueprint, 'r') as f:
            regions = [x.strip().split() for x in f.readlines()]
        env.log('{:,d} families with a total of {:,d} samples will be scanned for {:,d} pre-defined units'.\
                format(len(data.families), len(data.samples), len(regions)))
        env.jobs = max(min(args.jobs, len(regions)), 1)
        regions.extend([None] * env.jobs)
        queue = Queue()
        try:
            faulthandler.enable(file=open(env.tmp_log + '.SEGV', 'w'))
            for i in regions:
                queue.put(i)
            jobs = [EncoderWorker(
                queue, len(regions), deepcopy(data),
                RegionExtractor(args.vcf, chr_prefix = args.chr_prefix, allele_freq_info = args.freq),
                MarkerMaker(args.size),
                LinkageWriter(len(samples_not_vcf))
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
            env.log('{:,d} units (from {:,d} variants) processed; '\
                '{:,d} Mendelian inconsistencies and {:,d} recombination events handled\n'.\
                format(env.success_counter.value,
                       env.variants_counter.value, 
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
                env.error("{:,d} super markers failed to be generated due to haplotyping failures!".\
                          format(env.chperror_counter.value))
            if fatal_errors:
                env.error("{:,d} or more super markers failed to be generated due to runtime errors!".\
                          format(fatal_errors))
            env.log('Archiving to directory [{}]'.format(env.cache_dir))
            cache.write(pres = [env.output], exts = ['.tped', '.freq'], files = [env.outputfam],
                        otherfiles = [args.vcf, args.tfam, args.blueprint])
    env.jobs = args.jobs
    # STEP 2: write to PLINK or mega2 format
    tpeds = [os.path.join(env.cache_dir, item) for item in os.listdir(env.cache_dir) if item.startswith(env.output) and item.endswith('.tped')]
    for fmt in args.format:
        env.log('Saving data to directory [{}] ...'.format(fmt.upper()))
        format_linkage(tpeds, env.outputfam, fmt, args.prevalence, args.inherit_mode)
    #if 'plink' in args.format:
    #    env.log('Saving data to directory [PLINK] ...')
    #    formatPlink(tpeds, [env.outputfam] * len(tpeds), 'PLINK')
    #if 'mega2' in args.format:
    #    env.log('Saving data to directory [MEGA2] ...')
    #    formatMega2('MEGA2/{}*'.format(env.output))
    #if 'mlink' in args.format:
    #    env.log('Saving data to directory [MLINK] ...')
    #    format_linkage(tpeds, [env.outputfam] * len(tpeds), 'MLINK')
    #runMlink(args.blueprint)
    #plotMlink()
    #please implement this section
    # STEP final: clean up unwanted files
    #for item in args.format:
    #    removeFiles(item.upper(), exclude = env.formats[item])
    #removeFiles(env.cache_dir, exclude = ['.cache'])
