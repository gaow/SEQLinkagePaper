#!/usr/bin/env python
# Author: Biao Li (biaol@bcm.edu) and Gao Wang
# Date: 03-01-2014
# Purpose: simulation program to evaluate statistical power of SEQLinkage algorithm 

import argparse, random, tempfile, os, sys, shutil, time, glob, tarfile, copy
import progressbar
import collections, csv
import numpy as np
import glob

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

def getColumn(fn, num, delim = None, exclude = None):
    if num > 0:
        num = num - 1
    with open(fn) as inf:
        output = []
        for line in inf:
            parts = line.split(delim) if delim is not None else line.split()
            if len(parts) > num and parts[num] != exclude:
                output.append(parts[num])
    return output

OFF_PROP_2more = {2: 0.6314, 3: 0.2556, '4_more': 0.1130}
OFF_PROP_2and3 = {2: 0.7118, 3: 0.2882}
OFF_PROP_3more = {3:0.6934, '4_more':0.3066}
GENEMAP = {'11':'MYO7A', '13':'GJB2', '7':'SLC26A4', '22':'MYH9'}

def arguments(parser):
    parser.add_argument('-g', '--genes',
                        type=str,
                        metavar='FILE',
                        nargs=2,
                        help='''Specify input gene file names (gene1 gene2)
                        ''')
    parser.add_argument('-s', '--offspring',
                        type=int,
                        metavar='INT',
                        nargs=2,
                        default=[3,8],
                        help='''Specify a range of allowed number of offspring each parent may produce, allowed minimum number is between 2 and 4, maximum between 5 and 10. Default to "3 8"''')
    parser.add_argument('-m', '--mode',
                        type=str,
                        default='compound_recessive',
                        choices=['recessive', 'dominant', 'compound_recessive'],
                        help='''Specify mode of inheritance of disease locus. Compound recessive is the so-called compound heterozygosity.''')
    parser.add_argument('-n', '--num-families',
                        dest = "numfamilies",
                        type=int,
                        metavar='INT',
                        default=1,
                        help='''Specify number of families''')
    parser.add_argument('-a', '--allelic-heterogeneous',
                        action='store_true',
                        dest = 'allelichet',
                        help='''Whether or not the causal variant for a gene (within the same loci) is different for different families.''')
    parser.add_argument('-p', '--props-locus-heterogeneity',
                        dest = "locusheterogenprop",
                        type=float,
                        metavar='FLOAT',
                        nargs=2,
                        default=[0.5,0.5],
                        help='''Specify proportion of locus heterogeneity between two genes in the sample (the two proportion should sum to 1.0), e.g. 0.5 0.5''')
    parser.add_argument('--recrate',
                        type=float,
                        default=0.01,
                        help='''Recombination rate on the gene region''')
    parser.add_argument('-r', '--num-replicates',
                        dest = 'numreps',
                        type=int,
                        metavar='INT',
                        default=1,
                        help='''Specify desired number of replicates''')
    parser.add_argument('--ofile',
                        dest = 'outfile',
                        type=str,
                        default='simulation',
                        help='''Specify output file path, simulated data will be save to a *.fam file in linkage format and a *.vcf file in variant call format''')
    parser.add_argument('--seed',
                        type=float,
                        default=None,
                        help='''Specify seed for random number generator, if left unspecified the current system time will be used''')
    parser.add_argument('--save',
                        action='store_true',
                        default=False,
                        help='''Save ped and vcf files for each replicate''')
    parser.add_argument('--ld',
                        action='store_true',
                        default=False,
                        help='''Fake LD structure on one var site adjacent to the causal var site''')
    parser.add_argument('--alpha',
                        type=float,
                        nargs=2,
                        default = [3.3,3.6],
                        help='''Significance levels for LOD and HLOD (LOD HLOD), default to 3.3 3.6
                        ''')
    parser.add_argument('--debug',
                        action='store_true',
                        help='''Turn on Debug mode in which progress bar will not be shown and if an error occurs detailed error messages will be printed''')
    parser.add_argument('--sim-only',
                        action='store_true',
                        dest = 'sim_only',
                        help='''Simulation only''')
    parser.add_argument('--print-header',
                        action='store_true',
                        dest = 'print_header',
                        help='''Print output table header and quit.''')


def main(args, unknown_args):
    '''
    main func: pass cmd input and simulate pedigree samples
    '''
    def cleanup(items):
        for item in ['{}.{}'.format(args.outfile, x) for x in items]:
            os.system('/bin/rm -rf {}'.format(item))
    ## check user input
    checkInput(args)
    ## print header
    if args.print_header:
        print("gene1,gene2,plod1,nlod1,phlod1,nhlod1,plod2,nlod2,phlod2,nhlod2,prop1,prop2,moi,ahet,fam_size,min_offspring,max_offspring")
        sys.exit()
    ## set seed
    if not args.seed:
        args.seed = time.time()
    random.seed(args.seed)
    ## update proportions of number of offspring
    offNumProp = updateOffNumProp(args.offspring)
    ## parse input gene info
    gene1, gene2 = parseGeneInfo(args.genes[0]), parseGeneInfo(args.genes[1])
    ## clean up directory
    cleanup(['lods','hlods'])
    ## simulation
    counter = {} 
    alphas = {'lods':args.alpha[0], 'hlods':args.alpha[1]}
    if not args.debug:
        pbar = progressbar.ProgressBar(widgets=['Simulation [{}]'.format(args.seed), ' ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA(), ' '], maxval=int(args.numreps)).start()
    for i in xrange(1, args.numreps+1):
        samples = []
        if args.allelichet:
            diseaseVariantIndices = None
        else:
            diseaseVariantIndices = [gene1['d_idx'][weightedRandomIdx(gene1['cumuProbs_dMaf'])],
                                     gene2['d_idx'][weightedRandomIdx(gene2['cumuProbs_dMaf'])]]
        for j in xrange(args.numfamilies):
            numOffspring = getNumOffspring(offNumProp)
            pedInfo = simPedigree([gene1, gene2], numOffspring, args.mode, args.locusheterogenprop,
                                  diseaseVariantIndices, familyID=j+1, recRate=args.recrate, fakeLD=args.ld)
            samples.extend(pedInfo)
        # write *.fam file per sample for pedigree structure info only
        # write *.vcf file per sample for variant info
        dirName, baseName = os.path.dirname(args.outfile), os.path.basename(args.outfile)
        fam = os.path.join(dirName, baseName + 'rep'+str(i) + ".fam")
        vcf = os.path.join(dirName, baseName + 'rep'+str(i) + ".vcf")
        writeVCF(samples, writePedsToFile(samples, fam, pedStructOnly=True), gene1, gene2, vcf)
        # also write *.ped file if args.save is true
        if args.save:
            ped = os.path.join(dirName, baseName + 'rep'+str(i) + ".ped")
            writePedsToFile(samples, ped, pedStructOnly=False)
        if args.sim_only:
            if not args.debug:
                pbar.update(i)
            continue
        if not args.save:
            cleanup(['vcf.gz', 'vcf.gz.tbi'])
        vcf = indexVCF(vcf, verbose = False)
        # linkage analysis
        cmd = "seqlink --vcf {} --fam {} --output {} {} 2> /dev/null".\
                  format(vcf, fam, args.outfile, " ".join(unknown_args))
        os.system(cmd)
        res = {'lods':{}, 'hlods':{}}
        for score in ['lods', 'hlods']:
            for fn in glob.glob('{}/heatmap/*.{}'.format(args.outfile, score)):
                for marker, value in zip(getColumn(fn, 1), getColumn(fn, -1)):
                    value = float(value)
                    # convert single SNV marker to gene marker
                    if ":" in marker:
                        marker = GENEMAP[marker.split(":")[0]]
                    if marker not in res[score]:
                        res[score][marker] = value
                    else:
                        if value > res[score][marker]:
                            res[score][marker] = value
        # write result to file and calculate significance
        for score in res:
            with open(args.outfile + '.{}'.format(score), 'a') as f:
                for marker in res[score]:
                    f.write('{}\t{}\n'.format(marker, res[score][marker]))
                    if marker not in counter:
                        counter[marker] = {'lods': [0, 0], 'hlods': [0, 0]}
                    if res[score][marker] >= alphas[score]:
                        counter[marker][score][0] += 1 
                        counter[marker][score][1] += 1 
                    else:
                        counter[marker][score][1] += 1 
        if not args.debug:
            pbar.update(i)
    if not args.debug:
        pbar.finish()    
    # report power calculation
    gs = [x.split('.')[0] for x in args.genes]
    out = [gs[0], gs[1],
           counter[gs[0]]['lods'][0] / float(counter[gs[0]]['lods'][1] + 1E-10) if gs[0] in counter else -9,
           counter[gs[0]]['lods'][1] if gs[0] in counter else -9,
           counter[gs[0]]['hlods'][0] / float(counter[gs[0]]['hlods'][1] + 1E-10) if gs[0] in counter else -9,
           counter[gs[0]]['hlods'][1] if gs[0] in counter else -9,
           counter[gs[1]]['lods'][0] / float(counter[gs[1]]['lods'][1] + 1E-10) if gs[1] in counter else -9,
           counter[gs[1]]['lods'][1] if gs[1] in counter else -9,
           counter[gs[1]]['hlods'][0] / float(counter[gs[1]]['hlods'][1] + 1E-10) if gs[1] in counter else -9,
           counter[gs[1]]['hlods'][1] if gs[1] in counter else -9,
           args.locusheterogenprop[0], args.locusheterogenprop[1], args.mode,
           int(args.allelichet), args.numfamilies,
           args.offspring[0], args.offspring[1]]
    print ','.join(map(str, out))
    return 0
    

def writeVCF(samples, sample_names, gene1, gene2, fileName):
    '''
    write variant info to 'filename' in vcf format
    '''
    fi = open(fileName, 'w')
    varInfo = np.transpose(np.array(samples)[:,6:])
    chrInfo = gene1['chr'] + gene2['chr']
    posInfo = gene1['pos'] + gene2['pos']
    numVars = len(varInfo)/2
    fi.write("##fileformat=VCFv4.0\n")
    fi.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sample_names) + '\n')
    for idx in xrange(numVars):
        fi.write('\t'.join([str(chrInfo[idx]), str(posInfo[idx]), 'V'+str(idx+1), 'A', 'C', '.', 'PASS', '.', 'GT']) + '\t')
        fi.write('\t'.join(str(i)+'/'+str(j) for i,j in zip(varInfo[2*idx], varInfo[2*idx+1])) + '\n') 
    fi.close()
    return
    
    
def writePedsToFile(peds, fileName, pedStructOnly=False):
    '''
    write pedigree samples to 'fileName' in linkage format,
    if 'pedStructOnly' is True, output only the first 6 columns
    '''
    sample_names = []
    peds_clone = copy.deepcopy(peds)
    fi = open(fileName, 'w')
    for ind_info in peds_clone:
        for i in [1,2,3]:
            if ind_info[i]:
                ind_info[i] = "{}:{}".format(ind_info[0], ind_info[i])
        if pedStructOnly:
            fi.write(' '.join(str(x) for x in ind_info[:6]) + '\n')
        else:
            fi.write(' '.join(str(x) for x in ind_info) + '\n')
        sample_names.append(ind_info[1])
    fi.close()
    return sample_names 
    
    
def checkInput(args):
    '''
    validate user inputs and raise InputError in case of exception
    '''
    if args.print_header:
        return
    # minimum # offspring between [2,4], maximum # offspring between [3,10]
    if args.offspring[0] not in [2,3,4] or args.offspring[1] not in range(3,11) or args.offspring[0] > args.offspring[1]:
        sys.exit('Inappropriate input of argument "offspring"')
    if not args.genes:
        sys.exit("Please specify input files!")
    for item in args.genes:
        if not os.path.isfile(item):
            sys.exit('Cannot find file {}'.format(item)) 
        if os.path.basename(item).split('.')[0] not in GENEMAP.values():
            sys.exit('Gene {} is not supported!'.format(item.split('.')[0]))
    return


def updateOffNumProp(offspringRange):
    '''
    return an updated proportion of number of offspring, according to both OFF_PROP and args.offspring
    '''
    minimum, maximum = offspringRange[0], offspringRange[1]
    offProp, lenRange = collections.OrderedDict({}), maximum-minimum+1
    # special cases
    if lenRange == 1:
        return offProp.update({offspringRange[0]:1.0})
    elif offspringRange == [2,3]:
        offProp.update({2: 0.7118, 3: 0.2882})
    elif offspringRange == [2,4]:
        offProp.update({2: 0.6314, 3: 0.2556, 4: 0.1130})
    elif offspringRange == [3,4]:
        offProp.update({3: 0.6934, 4:0.3066})
    #
    else:
        if minimum == 2:
            Sn = 0.1130
            offProp.update({2: 0.6314, 3: 0.2556})
        elif minimum == 3:
            Sn = 0.3066
            offProp.update({3: 0.6934})
        else: # minimum == 4
            Sn = 1.0
        p = 2/3.*Sn/(1-(1./3**(maximum-4+1)))
        [offProp.update({n:p*((1/3.)**i)}) for i,n in enumerate(range(4, maximum+1))]
    return offProp
    


def parseGeneInfo(fileName):
    '''
    parse input gene file (*.tsv) and return dict obj of gene info
    '''
    info = {'pos':[], 'maf':[], 'annotation':[], 'd_idx':[], 'd_maf':[]}
    try:
        reader = csv.reader(open(fileName, 'rb'), delimiter='\t')
        rows = list(reader)
    except Exception:
        raise ValueError("Inappropriate argument value --genes")
    info['chr'] = [int(i[0]) for i in rows]
    info['pos'] = [int(i[1]) for i in rows]
    info['maf'] = [float(i[2]) for i in rows]
    info['annotation'] = [bool(int(i[4])) for i in rows]
    info['d_idx'] = [i for i,j in enumerate(info['annotation']) if j]
    info['nd_idx'] = [i for i,j in enumerate(info['annotation']) if not j]
    info['d_maf'] = [info['maf'][i] for i in info['d_idx']]
    cumu_dMaf = [sum(info['d_maf'][:i]) for i in range(1, len(info['d_maf'])+1)]
    sum_dMaf = sum(info['d_maf'])
    info['cumuProbs_dMaf'] = [i/sum_dMaf for i in cumu_dMaf]
    return info


def simPedigree(genes, numOffspring, mode, hetero, dVarIndices, familyID, recRate=0.01, fakeLD=False):
    '''
    simulate two-generational pedigree based on given input gene info, number of offspring, mode of inheritance and allelic heterogenity ratio.
    Note: at least two affected offspring are required per family sample
    '''
    causalGeneIdx = 0 if random.random() < hetero[0] else 1
    markerGeneIdx = 1 if causalGeneIdx == 0 else 0
    if dVarIndices is None:
        dVarIndices = [genes[0]['d_idx'][weightedRandomIdx(genes[0]['cumuProbs_dMaf'])],
                       genes[1]['d_idx'][weightedRandomIdx(genes[1]['cumuProbs_dMaf'])]]
    dVarIdx = dVarIndices[causalGeneIdx]
    #
    causalHaps, diseaseStatus = getCausalHaps(genes[causalGeneIdx], mode, numOffspring, dVarIdx, recRate, fakeLD)
    markerHaps = getMarkerHaps(genes[markerGeneIdx], numOffspring, recRate)
    pedInfo = createPedInfoDict(causalGeneIdx, causalHaps, markerHaps, diseaseStatus, familyID)
    return pedInfo


def createPedInfoDict(causalGeneIdx, causalHaps, markerHaps, diseaseStatus, familyID):
    pedInfo = []
    gene1Haps = causalHaps if causalGeneIdx == 0 else markerHaps
    gene2Haps = markerHaps if causalGeneIdx == 0 else causalHaps
    for idx, indID in enumerate(gene1Haps.keys()):
        tmp = [familyID, indID]
        if indID == 1:
            #sex1 = 1 if random.random() < 0.5 else 2
            tmp += [0,0,1]
        elif indID == 2:
            #sex2 = 1 if sex1 == 2 else 2
            tmp += [0,0,2]
        else:
            sex = 1 if random.random() < 0.5 else 2
            tmp += [1, 2, sex]
        tmp += [diseaseStatus[idx]]
        [tmp.extend([v1, v2]) for v1, v2 in zip(gene1Haps[indID][0], gene1Haps[indID][1])]
        [tmp.extend([v1, v2]) for v1, v2 in zip(gene2Haps[indID][0], gene2Haps[indID][1])]
        pedInfo.append(tmp)
    return pedInfo


def getMarkerHaps(geneInfo, numOffspring, recRate):
    '''
    simulate non-disease-causing haplotypes
    '''
    markerHaps = collections.OrderedDict({})
    parHaps = [genMarkerHap(geneInfo) for _ in range(4)]
    markerHaps[1], markerHaps[2] = parHaps[:2], parHaps[2:]
    if not recRate > 0: 
        [markerHaps.update({i+3:[parHaps[random.choice([0,1])], parHaps[random.choice([2,3])]]}) for i in range(numOffspring)]
    else:
        for idx in range(numOffspring):
            indHap1 = genGametHap(parHaps[0], parHaps[1], recRate)
            indHap2 = genGametHap(parHaps[2], parHaps[3], recRate)
            markerHaps[idx+3] = [indHap1, indHap2]
    return markerHaps
    

def getCausalHaps(geneInfo, mode, numOffspring, dVarIdx, recRate, fakeLD=False):
    '''
    Allowed dominant parental haplotype patterns:
    +-/--, ++/--, +-/+-
    Allowed recessive parental haplotype patterns:
    +-/+-, ++/+-
    '''
    causalHaps = collections.OrderedDict({})
    parentalHapCausality = [None,None,None,None] # 1 - causal; 0 - non-causal
    parAff = [None,None] # 1 - unaffected; 2 - affected
    avoids = [dVarIdx] if fakeLD else []
    while True:
        hap2, parentalHapCausality[1] = genCausalHap(geneInfo)
        hap4, parentalHapCausality[3] = genCausalHap(geneInfo)  
        if 'dominant' in mode:
            hap1, parentalHapCausality[0] = genCausalHap(geneInfo, addCausalVars=[dVarIdx])
            hap3, parentalHapCausality[2] = genCausalHap(geneInfo)
            # up to two '+' are allowed
            if parentalHapCausality.count(1) > 2:
                continue
            parAff[0] = 2
            parAff[1] = 2 if sum(parentalHapCausality[2:]) > 0 else 1
        elif 'recessive' in mode:
            hap1, parentalHapCausality[0] = genCausalHap(geneInfo, addCausalVars=[dVarIdx])
            if '_' in mode:
                # compound recessive. randomly choose another dVarIdx to use for this family instead of using the provided one
                causalVarIdx = geneInfo['d_idx'][weightedRandomIdx(geneInfo['cumuProbs_dMaf'], avoid = geneInfo['d_idx'].index(dVarIdx))]
                hap3, parentalHapCausality[2] = genCausalHap(geneInfo, addCausalVars=[causalVarIdx])
                if fakeLD:
                    avoids.append(causalVarIdx)
            else:
                hap3, parentalHapCausality[2] = genCausalHap(geneInfo, addCausalVars=[dVarIdx])
            # up to three '+' are allowed
            if parentalHapCausality.count(1) > 3:
                continue
            parAff[0] = 2 if sum(parentalHapCausality[:2]) == 2 else 1
            parAff[1] = 2 if sum(parentalHapCausality[2:]) == 2 else 1
        else:
            pass
        if not recRate > 0: # without recombination
            offspringAff, offHapIdx = genOffspringAffAndHapIdx(mode, parentalHapCausality,
                                                               numOffspring, [hap1, hap2, hap3, hap4],
                                                               geneInfo['d_idx'])
            if offspringAff.count(2) >= 2:
                # fill in parental and offspring causalHaps dict (1 - father; 2 - mother; 3,...,n - offspring)
                causalHaps[1], causalHaps[2] = [hap1, hap2], [hap3, hap4]
                parHaps = [hap1, hap2, hap3, hap4]
                [causalHaps.update({i+3:[parHaps[j[0]], parHaps[j[1]]]}) for i,j in enumerate(offHapIdx)]
                break
        else: # with recombination
            offspringAff, offspringHaps = genOffspringAffAndHaps(mode, numOffspring, [hap1, hap2, hap3, hap4], geneInfo['d_idx'], recRate)
            if offspringAff.count(2) >= 2:
                causalHaps[1], causalHaps[2] = [hap1, hap2], [hap3, hap4]
                [causalHaps.update({i+3:x}) for i,x in enumerate(offspringHaps)]
                break
    aff = parAff + offspringAff
    # fake LD structure on causal haps if fakeLD is True
    if fakeLD:
        listIdx = copy.deepcopy(geneInfo['d_idx'])
        if len(avoids) == 2:
            listIdx.remove(avoids[-1])
        tmpIdx = listIdx.index(dVarIdx)
        if tmpIdx == 0:
            fakeIdx = listIdx[1]
        elif tmpIdx == len(listIdx) - 1:
            fakeIdx = listIdx[-2]
        else:
            fakeIdx = listIdx[tmpIdx-1 if random.random() < 0.5 else tmpIdx+1]
        for (hap1, hap2) in causalHaps.values():
            hap1[fakeIdx] = hap1[dVarIdx]
            hap2[fakeIdx] = hap2[dVarIdx]
    #print dVarIdx
    #print parentalHapCausality
    #print offspringAff, offHapIdx
    #if parentalHapCausality.count(1) == 2:
    return causalHaps, aff
        

def genOffspringAffAndHaps(mode, numOffspring, haps, d_idx, recRate):
    def check_recessive_causal(h1, h2):
        if 2 in [x+y for i,(x,y) in enumerate(zip(h1, h2)) if i in d_idx]:
            return 2
        else:
            return 1
    #
    offHaps = [None] * numOffspring
    aff = [1] * numOffspring
    for idx in range(numOffspring):
        indHap1 = genGametHap(haps[0], haps[1], recRate)
        indHap2 = genGametHap(haps[2], haps[3], recRate)
        offHaps[idx] = [indHap1, indHap2]
        g = [checkIfCausal(d_idx, indHap1), checkIfCausal(d_idx, indHap2)]
        n = g.count(True)
        if 'dominant' in mode:
            aff[idx] = 2 if n > 0 else 1
        elif 'recessive' in mode:
            if '_' in mode:
                aff[idx] = 2 if n == 2 else 1
            else:
                aff[idx] = 1 if n < 2 else check_recessive_causal(indHap1, indHap2)
        else:
            raise ValueError("Inappropriate argument value '--mode'")
    return aff, offHaps


def genOffspringAffAndHapIdx(mode, parentalHapCausality, numOffspring, haps, d_idx):
    def check_recessive_causal(idx):
        if 2 in [x + y for i, (x, y) in enumerate(zip(haps[hapIdx[idx][0]], haps[hapIdx[idx][1]]))
                 if i in d_idx]:
            return 2
        else:
            return 1
    #
    hapIdx = [None] * numOffspring
    aff = [1] * numOffspring # 1 - unaffected; 2 - affected
    for idx in range(numOffspring):
        hapIdx[idx] = [random.choice([0,1]), random.choice([2,3])]
        g = [parentalHapCausality[hapIdx[idx][0]], parentalHapCausality[hapIdx[idx][1]]]
        n = g.count(1)
        if 'dominant' in mode:
            aff[idx] = 2 if n > 0 else 1
        elif 'recessive' in mode:
            if '_' in mode:
                aff[idx] = 2 if n == 2 else 1
            else:
                aff[idx] = 1 if n == 1 else check_recessive_causal(idx)
        else:
            raise ValueError("Inappropriate argument value '--mode'")
    return aff, hapIdx
            

def genGametHap(hap1, hap2, recRate):
    '''
    generate a gamet haplotype allowing for recombination
    '''
    if random.random() < recRate:
        spot = random.choice(range(len(hap1))[1:])
        hap3 = hap1[:spot] + hap2[spot:]
        hap4 = hap2[:spot] + hap1[spot:]
        return random.choice([hap1, hap2, hap3, hap4])
    else:
        return random.choice([hap1, hap2])


def genMarkerHap(geneInfo):
    hap = [0] * len(geneInfo['maf'])
    for idx in geneInfo['nd_idx']:
        if random.random() < geneInfo['maf'][idx]:
            hap[idx] = 1
    return hap
  

def genCausalHap(geneInfo, addCausalVars=[]):
    '''
    generate 1st parental haplotype that contains 1 disease variant
    '''
    ifCausal = False
    hap = [0] * len(geneInfo['maf'])
    for idx in xrange(len(geneInfo['maf'])):
        if random.random() < geneInfo['maf'][idx]:
            hap[idx] = 1
    # add causal variants on given sites
    for idx in addCausalVars:
        hap[idx] = 1
    if len(addCausalVars) > 0 or checkIfCausal(geneInfo['d_idx'], hap):
        ifCausal = True
    return hap, ifCausal


def checkIfCausal(causalSites, hap):
    for d_idx in causalSites:
        if hap[d_idx] == 1:
            return True
    return False


def getNumOffspring(offNumProp):
    '''
    return randomly generated number of offspring according to proportion of number of offspring
    '''
    nums, probs = offNumProp.keys(), offNumProp.values()
    cumuProbs = [sum(probs[:i]) for i in range(1, len(probs)+1)]
    idx = weightedRandomIdx(cumuProbs)
    return nums[idx]


def weightedRandomIdx(cumuProbs, avoid = None):
    '''
    return a weighted random choice 
    '''
    if len(cumuProbs) == 1:
        avoid = None
    try:
        len(avoid)
    except:
        avoid = [avoid]
    while True:
        randNum = random.random()
        for idx, p in enumerate(cumuProbs):
            if randNum < p and idx not in avoid:
                return idx


if __name__ == '__main__':
    master_parser = argparse.ArgumentParser(
        description = '''Program to generate two generational family samples with two gene regions and perform power calculation with SEQLinkage program; All unknown args will be passed to seqlink program''',
        prog = 'seqlink-pc',
        epilog = '''Biao Li (biaol@bcm.edu) and Gao Wang (c) 2014.'''
    )
    master_parser.add_argument('--version,', action='version', version='%(prog)s 0.1.0')
    arguments(master_parser)
    master_parser.set_defaults(func=main)
    # getting arguments
    args, unknown_args = master_parser.parse_known_args()
    #print vars(args)
    # calling associated functions
    if args.debug:
        args.func(args, unknown_args)
    else:
        try:
            args.func(args, unknown_args)
        except Exception as e:
            sys.exit("An ERROR has occured: {}".format(e))
