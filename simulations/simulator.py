#!/usr/bin/env python
# Author: Biao Li (biaol@bcm.edu)
# Date: 03-01-2014
# Purpose: simulation program to test SEQLinco 


import argparse, random, tempfile, os, sys, shutil, time, glob, tarfile
import progressbar
import collections, csv

OFF_PROP_2more = {2: 0.6314, 3: 0.2556, '4_more': 0.1130}
OFF_PROP_2and3 = {2: 0.7118, 3: 0.2882}
OFF_PROP_3more = {3:0.6934, '4_more':0.3066}

#GENE_1 = "/Users/biao/Projects/SEQLinco/simulations/GJB2.tsv"
#GENE_2 = "/Users/biao/Projects/SEQLinco/simulations/SLC26A4.tsv"
GENE_1 = "GJB2.tsv"
GENE_2 = "SLC26A4.tsv"

def arguments(parser):
    parser.add_argument('-g', '--genes',
                        type=str,
                        metavar='FILE, FILE',
                        nargs='+',
                        default=[GENE_1, GENE_2],
                        help='''Specify input gene file names (gene1 gene2)
                        ''')
    parser.add_argument('-o', '--offspring',
                        type=int,
                        metavar='INT, INT',
                        nargs='+',
                        default=[2,5],
                        help='''Specify a range of allowed number of offspring each parent may produce, e.g. 2 5, allowed minimum number is between 2 and 4, maximum between 5 and 10''')
    parser.add_argument('-m', '--mode',
                        type=str,
                        default='dominant',
                        choices=['recessive', 'dominant', 'compound_recessive', 'compound_dominant'],
                        help='''Specify mode of inheritance of disease locus, e.g. recessive or dominant or compound_recessive or compound_dominant''')
    parser.add_argument('-s', '--samplesize',
                        type=int,
                        metavar='INT',
                        default=2,
                        help='''Specify sample size, number of families''')
    parser.add_argument('-a', '--allelicheteroprop',
                        type=float,
                        metavar='FLOAT, FLOAT',
                        nargs='+',
                        default=[0.5,0.5],
                        help='''Specify proportion of allelic heterogeneity between two genes in the sample, e.g. 0.5 0.5''')
    parser.add_argument('-r', '--numreps',
                        type=int,
                        metavar='INT',
                        default=3,
                        help='''Specify desired number of replicates (sample size)''')
    parser.add_argument('-f', '--outfile',
                        type=str,
                        default='sample',
                        help='''Specify output file name, simulated data will be save to a *.ped file in linkage format and a *.vcf file in variant seq format''')
    parser.add_argument('--seed',
                        type=float,
                        default=None,
                        help='''Specify seed for random number generator, if left unspecified the current system time will be used''')
    parser.add_argument('-d', '--debug',
                        type=bool,
                        default=False,
                        help='''Turn on Debug mode in which progress bar will not be shown and if an error occurs detailed error messages will be printed''')


def main(args):
    '''
    main func: pass cmd input and simulate pedigree samples
    '''
    ## check user input
    checkInput(args)
    ## set seed
    if args.seed:
        random.seed(args.seed)
    else:
        random.seed(time.time())
    ## update proportions of number of offspring
    offNumProp = updateOffNumProp(args.offspring)
    ## parse input gene info
    gene1, gene2 = parseGeneInfo(args.genes[0]), parseGeneInfo(args.genes[1])
    ## simulation
    if not args.debug:
        pbar = progressbar.ProgressBar(widgets=['Simulating for {} replicates'.format(args.numreps), ' ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA(), ' '], maxval=int(args.numreps)).start()
    for i in xrange(1, args.numreps+1):
        # per replicate
        samples = []
        diseaseVariantIndices = [weightedRandomIdx(gene1['cumuProbs_dMaf']), weightedRandomIdx(gene2['cumuProbs_dMaf'])]
        for j in xrange(args.samplesize):
            numOffspring = getNumOffspring(offNumProp)
            pedInfo = simPedigree([gene1, gene2], numOffspring, args.mode, args.allelicheteroprop, diseaseVariantIndices, familyID=j+1)
            samples.extend(pedInfo)
        # write *.ped file per sample for pedigree structure info only
        writePedsToFile(samples, args.outfile+"_rep"+str(i)+".ped", pedStructOnly=True)
        # write *.vcf file per sample for variant info
        writeVCF(samples, args.outfile+"_rep"+str(i)+".vcf")
        if not args.debug:
            pbar.update(i)
    if not args.debug:
        pbar.finish()    
    #
        # FIXME! here do something about simulated pedigree samples,
        #### FIXME!
        # e.g. save to file, do analysis, etc...
        # saveToFile(...)
    return
    

def writeVCF(samples, fileName):
    pass
    
    
def writePedsToFile(peds, fileName, pedStructOnly=False):
    '''
    write pedigree samples to 'fileName' in linkage format,
    if 'pedStructOnly' is True, output only the first 6 columns
    '''
    if not fileName.endswith('.ped'):
        fileName += '.ped'
    fi = open(fileName, 'w')
    for ind_info in peds:
        if pedStructOnly:
            fi.write(' '.join(str(x) for x in ind_info[:6]) + '\n')
        else:
            fi.write(' '.join(str(x) for x in ind_info) + '\n')
    fi.close()
    return    
    
    
def checkInput(args):
    '''
    validate user inputs and raise InputError in case of exception
    '''
    # minimum # offspring between [2,4], maximum # offspring between [3,10]
    if args.offspring[0] not in [2,3,4] or args.offspring[1] not in range(3,11) or args.offspring[0] > args.offspring[1]:
        raise ValueError('Inappropriate input of argument "offspring"')
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


def simPedigree(genes, numOffspring, mode, hetero, dVarIndices, familyID):
    '''
    simulate two-generational pedigree based on given input gene info, number of offspring, mode of inheritance and allelic heterogenity ratio.
    Note: at least two affected offspring are required per family sample
    '''
    causalGeneIdx = 0 if random.random() < hetero[0] else 1
    markerGeneIdx = 1 if causalGeneIdx == 0 else 0
    dVarIdx = dVarIndices[causalGeneIdx]
    #
    causalHaps, diseaseStatus = getCausalHaps(genes[causalGeneIdx], mode, numOffspring, dVarIdx)
    markerHaps = getMarkerHaps(genes[markerGeneIdx], numOffspring)
    pedInfo = createPedInfoDict(causalGeneIdx, causalHaps, markerHaps, diseaseStatus, familyID)
    return pedInfo


def createPedInfoDict(causalGeneIdx, causalHaps, markerHaps, diseaseStatus, familyID):
    pedInfo = []
    gene1Haps = causalHaps if causalGeneIdx == 0 else markerHaps
    gene2Haps = markerHaps if causalGeneIdx == 0 else causalHaps
    for idx, indID in enumerate(gene1Haps.keys()):
        tmp = [familyID, indID]
        if indID == 1:
            sex1 = 1 if random.random() < 0.5 else 2
            tmp += [0,0,sex1]
        elif indID == 2:
            sex2 = 1 if sex1 == 2 else 2
            tmp += [0,0,sex2]
        else:
            sex = 1 if random.random() < 0.5 else 2
            tmp += [sex1,sex2,sex]
        tmp += [diseaseStatus[idx]]
        [tmp.extend([v1, v2]) for v1, v2 in zip(gene1Haps[indID][0], gene1Haps[indID][1])]
        [tmp.extend([v1, v2]) for v1, v2 in zip(gene2Haps[indID][0], gene2Haps[indID][1])]
        pedInfo.append(tmp)
    return pedInfo


def getMarkerHaps(geneInfo, numOffspring):
    '''
    simulate non-disease-causing haplotypes
    '''
    markerHaps = collections.OrderedDict({})
    parHaps = [genMarkerHap(geneInfo) for _ in range(4)]
    markerHaps[1], markerHaps[2] = parHaps[:2], parHaps[2:]
    [markerHaps.update({i+3:[parHaps[random.choice([0,1])], parHaps[random.choice([2,3])]]}) for i in range(numOffspring)]
    return markerHaps

def getCausalHaps(geneInfo, mode, numOffspring, dVarIdx):
    '''
    if mode in ['dominant', 'recessive'], use given dVarIdx as disease variant indices for simulating the first causal haplotype, otherwise re-generate such indices if mode in ['compound_dominant', 'compound_recessive']
    Allowed dominant parental haplotype patterns:
    +-/--, ++/--, +-/+-
    Allowed recessive parental haplotype patterns:
    +-/+-, ++/+-
    '''
    causalHaps = collections.OrderedDict({})
    parentalHapCausality = [1,0,0,0] # 1 - causal; 0 - non-causal
    parAff = [1,1] # 1 - unaffected; 2 - affected
    if mode in ['compound_dominant', 'compound_recessive']:
        # randomly choose another dVarIdx to use per family instead of using the provided one
        dVarIdx = weightedRandomIdx(geneInfo['cumuProbs_dMaf'])
    while True:
        hap2, parentalHapCausality[1] = genCausalHap(geneInfo)
        hap4, parentalHapCausality[3] = genCausalHap(geneInfo)  
        if mode in ['dominant', 'compound_dominant']:
            hap1, parentalHapCausality[0] = genCausalHap(geneInfo, addCausalVars=[dVarIdx])
            hap3, parentalHapCausality[2] = genCausalHap(geneInfo)
            # up to two '+' are allowed
            if parentalHapCausality.count(1) > 2:
                continue
            parAff[0] = 2
            parAff[1] = 2 if sum(parentalHapCausality[2:]) > 0 else 1
        elif mode in ['recessive', 'compound_recessive']:
            hap1, parentalHapCausality[0] = genCausalHap(geneInfo, addCausalVars=[dVarIdx])
            hap3, parentalHapCausality[2] = genCausalHap(geneInfo, addCausalVars=[dVarIdx])
            # up to three '+' are allowed
            if parentalHapCausality.count(1) > 3:
                continue
            parAff[0] = 2 if sum(parentalHapCausality[:2]) == 2 else 1
            parAff[1] = 2 if sum(parentalHapCausality[2:]) == 2 else 1
        else:
            pass
        offspringAff, offHapIdx = genOffspringAffAndHapIdx(mode, parentalHapCausality, numOffspring)
        if offspringAff.count(2) >= 2:
            # fill in parental and offspring causalHaps dict (1 - father; 2 - mother; 3,...,n - offspring)
            causalHaps[1], causalHaps[2] = [hap1, hap2], [hap3, hap4]
            parHaps = [hap1, hap2, hap3, hap4]
            [causalHaps.update({i+3:[parHaps[j[0]], parHaps[j[1]]]}) for i,j in enumerate(offHapIdx)]
            break
    aff = parAff + offspringAff
    #print dVarIdx
    #print parentalHapCausality
    #print offspringAff, offHapIdx
    #if parentalHapCausality.count(1) == 2:
    return causalHaps, aff
        

def genOffspringAffAndHapIdx(mode, parentalHapCausality, numOffspring):
    hapIdx = [None] * numOffspring
    aff = [1] * numOffspring # 1 - unaffected; 2 - affected
    for idx in range(numOffspring):
        hapIdx[idx] = [random.choice([0,1]), random.choice([2,3])]
        g = [parentalHapCausality[hapIdx[idx][0]], parentalHapCausality[hapIdx[idx][1]]]
        n = g.count(1)
        if mode in ['dominant', 'compound_dominant']:
            aff[idx] = 2 if n > 0 else 1
        elif mode in ['recessive', 'compound_recessive']:
            aff[idx] = 2 if n == 2 else 1
        else:
            raise ValueError("Inappropriate argument value '--mode'")
    return aff, hapIdx
            


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


def weightedRandomIdx(cumuProbs):
    '''
    return a weighted random choice 
    '''
    randNum = random.random()
    for idx, p in enumerate(cumuProbs):
        if randNum < p:
            return idx
    


if __name__ == '__main__':
    master_parser = argparse.ArgumentParser(
        description = '''Program to generate two generational family samples with two gene regions''',
        prog = 'simSEQLinco',
        epilog = '''Biao Li (biaol@bcm.edu) and Hang Dai (hang.dai@bcm.edu) (c) 2014.'''
    )
    master_parser.add_argument('--version,', action='version', version='%(prog)s 0.1.0')
    arguments(master_parser)
    master_parser.set_defaults(func=main)
    # getting arguments
    args = master_parser.parse_args()
    #print vars(args)
    # calling associated functions
    if args.debug:
        args.func(args)
    else:
        try:
            args.func(args)
        except Exception as e:
            sys.exit("An ERROR has occured: {}".format(e))

