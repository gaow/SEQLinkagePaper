#!/usr/bin/python
# Copyright (c) 2013, Gao Wang <gaow@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
from __future__ import print_function
from SEQLinkage import HOMEPAGE
from SEQLinkage.Utils import *
from SEQLinkage.Runner import *
from multiprocessing import Process, Queue
from collections import Counter, OrderedDict, defaultdict
import itertools
from zipfile import ZipFile
from copy import deepcopy
import sys, faulthandler

if sys.version_info.major == 2:
    from cstatgen import cstatgen_py2 as cstatgen
else:
    from cstatgen import cstatgen_py3 as cstatgen
from cstatgen.egglib import Align

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
    1. samples in FAM but not in VCF --> ERROR
    2. samples in VCF but not in FAM --> give a message'''
    a_not_b = list(set(samp1).difference(set(samp2)))
    b_not_a = list(set(samp2).difference(set(samp1)))
    if b_not_a:
        env.log('{:,d} samples found in FAM file but not in VCF file:\n{}'.\
                           format(len(b_not_a), '\n'.join(b_not_a)))
    if a_not_b:
        env.log('{:,d} samples in VCF file will be ignored due to absence in FAM file'.format(len(a_not_b)))
    return a_not_b, b_not_a

class Cache:
    def __init__(self, cache_dir, cache_name, params):
        self.cache_dir = cache_dir
        self.cache_name = os.path.join(cache_dir, cache_name + '.cache')
        self.cache_info = os.path.join(cache_dir, '.info.' + cache_name)
        self.param_info = os.path.join(cache_dir, '.conf.' + cache_name)
        mkpath(cache_dir)
        self.id = '.' 
        self.params = params
        self.infofiles = [params['vcf'], params['tfam'], params['blueprint']] if params['blueprint'] else [params['vcf'], params['tfam']]
        self.infofiles.append(self.cache_name)
        self.pchecklist = {'.vcf': ['bin', 'single_markers'],
                           '.linkage': ['prevalence', 'inherit_mode', 'wild_pen',
                                        'muta_pen', 'theta_max', 'theta_inc'],
                           '.analysis': ['prevalence', 'inherit_mode', 'wild_pen',
                                        'muta_pen', 'theta_max', 'theta_inc']}

    def setID(self, ID):
        self.id = "." + str(ID)
        
    def check(self):
        if not os.path.isfile(self.cache_info) or not os.path.isfile(self.param_info + self.id):
            return False
        with open(self.cache_info, 'r') as f:
            lines = [item.strip().split() for item in f.readlines()]
        for line in lines:
            if not os.path.isfile(line[0]) or line[1] != calculateFileMD5(line[0]):
                return False
        params = {}
        with open(self.param_info + self.id, 'r') as f:
            lines = [item.strip().split() for item in f.readlines()]
        for line in lines:
            params[line[0]] = line[1]
        for item in self.pchecklist[self.id]:
            if params[item] != str(self.params[item]):
                return False
        return True

    def load(self, target_dir = None, names = None):
        if target_dir is None:
            target_dir = self.cache_dir
        with ZipFile(self.cache_name) as f:
            if names is None:
                f.extractall(target_dir)
            else:
                for item in f.namelist():
                    if any([item.startswith(x) for x in names]):
                        f.extract(item, target_dir)

    def write(self, source_dir = None, arcroot = '/', pres = [], exts = [],
              files = [], mode = 'w'):
        '''Add files to cache'''
        if source_dir is None:
            source_dir = self.cache_dir
        with ZipFile(self.cache_name, mode) as f:
            if source_dir != self.cache_dir:
                zipdir(source_dir, f, arcroot = arcroot)
            else:
                for item in os.listdir(source_dir):
                    if ((any([item.endswith(x) for x in exts]) or len(exts) == 0) \
                        and (any([item.startswith(x) for x in pres]) or len(pres) == 0)) \
                        or item in files:
                        f.write(os.path.join(source_dir, item), arcname=os.path.join(arcroot, item))
        signatures = ['{}\t{}'.format(x, calculateFileMD5(x)) for x in self.infofiles if os.path.isfile(x)]
        with open(self.cache_info, 'w') as f:
            f.write('\n'.join(signatures))
        with open(self.param_info + self.id, 'w') as f:
            f.write('\n'.join(["{}\t{}".format(item, self.params[item]) for item in self.pchecklist[self.id]]))

    def clear(self, pres = [], exts = []):
        for fl in glob.glob(self.cache_info + "*") + [self.cache_name]:
            try:
                os.remove(fl)
            except OSError:
                pass
        #
        for pre, ext in itertools.product(pres, exts): 
            for fl in glob.glob(os.path.join(self.cache_dir, pre) +  "*" + ext): 
                try:
                    os.remove(fl)
                except OSError:
                    pass
        

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
        '''member is one line of FAM file, [fid, sid, pid, mid, sex, pheno]'''
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
                names.append("V{}-{}".format(idx, self.variants[idx][1]))
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
        self.chrom = self.startpos = self.endpos = self.name = None
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
        self.chrom, self.startpos, self.endpos, self.name = region[:4]
        self.startpos = int(self.startpos)
        self.endpos = int(self.endpos) + 1
        if self.chrom in ['X','23']:
            if self.xchecker.check(self.startpos) or self.xchecker.check(self.endpos): 
                self.chrom = 'XY'
        if self.chrom in ['Y','24']:
            if self.ychecker.check(self.startpos) or self.ychecker.check(self.endpos):
                self.chrom = 'XY'
        if self.chr_prefix and not self.chrom.startswith(self.chr_prefix):
            self.chrom = self.chr_prefix + self.chrom

    
class MarkerMaker:
    def __init__(self, wsize, maf_cutoff = None):
        self.missings = ("0", "0")
        self.gtconv = {'1':0, '2':1}
        self.haplotyper = cstatgen.HaplotypingEngine(verbose = env.debug)
        if wsize == 0 or wsize >= 1:
            self.r2 = None
        else:
            self.r2 = wsize
        self.coder = cstatgen.HaplotypeCoder(wsize)
        self.maf_cutoff = maf_cutoff

    def apply(self, data):
        # temp raw haplotype, maf and variant names data
        haplotypes = OrderedDict() 
        mafs = {}
        varnames = {}
        try:
            # haplotyping plus collect found allele counts
            # and computer founder MAFS
            self.__Haplotype(data, haplotypes, mafs, varnames)
            if len(varnames):
                if not any ([len(varnames[x]) - 1 for x in varnames]):
                    # all families have only one variant
                    self.__AssignSNVHaplotypes(data, haplotypes, mafs, varnames)
                else:
                    # calculate LD clusters using founder haplotypes
                    clusters = self.__ClusterByLD(data, haplotypes, varnames)
                    # recoding the genotype of the region
                    self.__CodeHaplotypes(data, haplotypes, mafs, varnames, clusters)
        except Exception as e:
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
        self.markers = ["V{}-{}".format(idx, item[1]) for idx, item in enumerate(data.variants)]
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
        #
        # Drop some variants if maf is greater than given threshold 
        #
        if self.maf_cutoff is not None:
            exclude_vars = []
            for v in mafs.keys():
                if mafs[v] > self.maf_cutoff:
                    exclude_vars.append(v)
            for i in haplotypes.keys():
                haplotypes[i] = listit(haplotypes[i])
                for j in range(len(haplotypes[i])):
                    haplotypes[i][j] = haplotypes[i][j][:2] + \
                      [x for idx, x in enumerate(haplotypes[i][j][2:]) if varnames[i][idx] not in exclude_vars]
                varnames[i] = [x for x in varnames[i] if x not in exclude_vars]
                # handle trivial data
                if len(varnames[i]) == 0:
                    for person in data.families[i]:
                        data[person] = self.missings
                    del varnames[i]
                    del haplotypes[i]
            # count how many variants are removed
            with env.commonvar_counter.get_lock():
                env.commonvar_counter.value += len(exclude_vars)


    def __ClusterByLD(self, data, haplotypes, varnames):
        if self.r2 is None:
            return None
        # get founder haplotypes
        founder_haplotypes = []
        markers = sorted(set(itertools.chain(*varnames.values())), key = lambda x: int(x.split("-")[0][1:]))
        for item in haplotypes:
            for ihap, hap in enumerate(haplotypes[item]):
                if not data.tfam.is_founder(hap[1]):
                    continue
                gt = [hap[2 + varnames[item].index(v)] if v in varnames[item] else '?' for v in markers]
                founder_haplotypes.append(("{}-{}".format(hap[1], ihap % 2), "".join([x[1] if x[0].isupper() else x[0] for x in gt])))
        # calculate LD blocks, use r2 measure
        ld = Align.create(founder_haplotypes).matrixLD(validCharacters="12")["r2"]
        blocks = []
        for j in ld:
            block = [j]
            for k in ld[j]:
                if ld[j][k] > self.r2:
                    block.append(k)
            if len(block) > 1:
                blocks.append(block)
        # get LD clusters
        return [[markers[idx] for idx in item] for item in list(connected_components(blocks))]
        

    def __CodeHaplotypes(self, data, haplotypes, mafs, varnames, clusters):
        # apply CHP coding
        if clusters is not None:
            clusters_idx = [[[varnames[item].index(x) for x in y] for y in clusters] for item in haplotypes]
        else:
            clusters_idx = [[[]] for item in haplotypes]
        self.coder.Execute(haplotypes.values(), [[mafs[v] for v in varnames[item]] for item in haplotypes], clusters_idx)
        if env.debug:
            with env.lock:
                self.coder.Print()
        # line: [fid, sid, hap1, hap2]
        for line in self.coder.GetHaplotypes():
            if not line[1] in data:
                # this sample is not in VCF file. Every variant site should be missing
                # they have to be skipped for now
                continue 
            data[line[1]] = (line[2].split(','), line[3].split(','))
            if len(data[line[1]][0]) > data.superMarkerCount:
                data.superMarkerCount = len(data[line[1]][0])
        # get MAF
        for item in haplotypes:
            data.maf[item] = self.coder.GetAlleleFrequencies(item)
        if env.debug:
            with env.lock:
                print("marker freqs = ", data.maf, "\n", file = sys.stderr)

                
    def __AssignSNVHaplotypes(self, data, haplotypes, mafs, varnames):
        for item in haplotypes:
            # each person's haplotype
            token = ''
            for idx, line in enumerate(haplotypes[item]):
                if not idx % 2:
                    token = line[2][1] if line[2][0].isupper() else line[2][0]
                else:
                    data[line[1]] = (token, line[2][1] if line[2][0].isupper() else line[2][0]) 
            # get maf
            data.maf[item] = [(1 - mafs[varnames[item][0]], mafs[varnames[item][0]])]
    

    def __FormatHaplotypes(self, data):
        # Reformat sample genotypes 
        for person in data:
            if type(data[person]) is not tuple:
                data[person] = self.missings
                continue
            diff = data.superMarkerCount - len(data[person][0])
            data[person] = zip(*data[person])
            if diff > 0:
                data[person].extend([self.missings] * diff)
        

class LinkageWriter:
    def __init__(self, num_missing_append = 0):
        self.chrom = self.prev_chrom = self.name = self.distance = self.distance_avg = self.distance_m = self.distance_f = None
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
            # genotypes
            gs = [data[s][0] for s in data.samples]
            if len(set(gs)) == 1:
                # everyone's genotype is the same (most likely missing or monomorphic)
                return 2
            self.tped += env.delimiter.join([self.chrom, self.name, self.distance, position] + \
                list(itertools.chain(*gs)) + self.missings*self.num_missing) + "\n"
            # freqs
            for k in data.maf:
                self.freq += env.delimiter.join([k, self.name] + map(str, data.maf[k][0])) + "\n" 
        else:
            # have to expand each region into mutiple chunks to account for different recomb points
            gs = zip(*[data[s] for s in data.samples])
            # sub-chunk id
            cid = 0
            skipped_chunk = []
            for idx, g in enumerate(gs):
                if len(set(g)) == 1:
                    skipped_chunk.append(idx)
                    continue
                cid += 1
                self.tped += \
                  env.delimiter.join([self.chrom, '{}[{}]'.format(self.name, cid), self.distance, position] + \
                  list(itertools.chain(*g)) + self.missings*self.num_missing) + "\n"
            if cid == 0:
                # everyone's genotype is the same (most likely missing or monomorphic)
                return 2
            # freqs
            for k in data.maf:
                cid = 0
                for idx in range(data.superMarkerCount):
                    if idx in skipped_chunk:
                        continue
                    if idx >= len(data.maf[k]):
                        break
                    cid += 1
                    self.freq += env.delimiter.join([k, '{}[{}]'.format(self.name, cid)] + \
                                                    map(str, data.maf[k][idx])) + "\n" 
        if self.counter < env.batch:
            self.counter += data.superMarkerCount
        else:
            self.commit()
        return 0

    def commit(self):
        if self.tped:
            with env.lock:
                with open(os.path.join(env.tmp_cache, '{}.chr{}.tped'.format(env.output, self.prev_chrom)),
                          'a') as f:
                    f.write(self.tped)
        if self.freq:
            with env.lock:
                with open(os.path.join(env.tmp_cache, '{}.chr{}.freq'.format(env.output, self.prev_chrom)),
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
        self.name, self.distance_avg, self.distance_m, self.distance_f = region[3:]
        self.distance = ";".join([self.distance_avg, self.distance_m, self.distance_f])


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
    downloadResources([('{}/uploads/genemap.txt'.format(HOMEPAGE), env.resource_dir),
                       ('{}/uploads/mlink'.format(HOMEPAGE), env.resource_bin),
                       ('{}/uploads/unknown'.format(HOMEPAGE), env.resource_bin),
                       ('{}/uploads/makeped'.format(HOMEPAGE), env.resource_bin),
                       ('{}/uploads/pedcheck'.format(HOMEPAGE), env.resource_bin)])
    cache = Cache(env.cache_dir, env.output, vars(args))
    cache.setID('vcf')
    # STEP 1: write encoded data to TPED format
    if not args.vanilla and cache.check():
        env.log('Loading cache data from archive ...')
        cache.load(target_dir = env.tmp_dir, names = ['CACHE'])
    else:
        # load VCF file header
        checkVCFBundle(args.vcf)
        cache.clear()
        try:
            vs = cstatgen.VCFstream(args.vcf)
        except Exception as e:
            env.error("{}".format(e), exit = True)
        samples_vcf = vs.GetSampleNames()
        if len(samples_vcf) == 0:
            env.error("Fail to extract samples from [{}]".format(args.vcf), exit = True)
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
        rewriteFamfile(os.path.join(env.tmp_cache, '{}.tfam'.format(env.output)),
                       data.tfam.samples, data.samples.keys() + samples_not_vcf)
        if args.single_markers:
            regions = [(x[0], x[1], x[1], "{}:{}".format(x[0], x[1]), '.', '.', '.')
                       for x in vs.GetGenomeCoordinates()]
            args.blueprint = None
        else:
            # load blueprint
            try:
                with open(args.blueprint, 'r') as f:
                    regions = [x.strip().split() for x in f.readlines()]
            except IOError:
                env.error("Cannot load regional marker blueprint [{}]. ".format(args.blueprint), exit = True)
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
                MarkerMaker(args.bin, maf_cutoff = args.maf_cutoff),
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
            if env.commonvar_counter.value:
                env.log('{:,d} variants ignored due to having MAF > {}'.\
                        format(env.commonvar_counter.value, args.maf_cutoff))
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
            cache.write(arcroot = 'CACHE', source_dir = env.tmp_cache)
    env.jobs = args.jobs
    # STEP 2: write to PLINK or mega2 format
    tpeds = [os.path.join(env.tmp_cache, item) for item in os.listdir(env.tmp_cache) if item.startswith(env.output) and item.endswith('.tped')]
    for fmt in args.format:
        cache.setID(fmt)
        if not args.vanilla and cache.check():
            env.log('Loading {} data from archive ...'.format(fmt.upper()))
            cache.load(target_dir = env.tmp_dir, names = [fmt.upper()])
        else:
            env.log('Generating {} format output ...'.format(fmt.upper()))
            format(tpeds, os.path.join(env.tmp_cache, "{}.tfam".format(env.output)),
                   args.prevalence, args.wild_pen, args.muta_pen, fmt,
                   args.inherit_mode, args.theta_max, args.theta_inc)
            cache.write(arcroot = fmt.upper(),
                        source_dir = os.path.join(env.tmp_dir, fmt.upper()), mode = 'a')
    mkpath(env.output)
    if args.run_linkage:
        cache.setID('analysis')
        if not args.vanilla and cache.check():
            env.log('Loading linkage analysis result from archive ...'.format(fmt.upper()))
            cache.load(target_dir = env.output, names = ['heatmap'])
        else:
            env.log('Running LINKAGE now ...')
            run_linkage(args.blueprint, args.theta_inc, args.theta_max, args.output_limit)
            cache.write(arcroot = 'heatmap', source_dir = os.path.join(env.output, 'heatmap'), mode = 'a') 
        html(args.theta_inc, args.theta_max, args.output_limit)
    else:
        env.log('Saving data to [{}] ...'.format(os.path.abspath(env.output)))
        cache.load(target_dir = env.output, names = [fmt.upper() for fmt in args.format])
