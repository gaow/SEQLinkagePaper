#!/usr/bin/env python
# Author: Biao Li (biaol@bcm.edu)
# Date: 03-01-2014
# Purpose: simulation program to test SEQLinco 


import argparse, random, tempfile, os, sys, shutil, time, glob, tarfile
import progressbar
import collections, csv

DEBUG = True

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
                        default='recessive',
                        choices=['recessive', 'dominant', 'compound_recessive', 'compound_dominant'],
                        help='''Specify mode of inheritance of disease locus, e.g. recessive or dominant or compound_recessive or compound_dominant''')
    parser.add_argument('-s', '--samplesize',
                        type=int,
                        metavar='INT',
                        default=3,
                        help='''Specify sample size, number of families''')
    parser.add_argument('-a', '--allelicheteroprop',
                        type=float,
                        metavar='FLOAT, FLOAT',
                        nargs='+',
                        default=[0.5,0.5],
                        help='''Specify proportion of allelic heterogeneity between two genes in the sample, e.g. 0.5 0.5''')
    parser.add_argument('-n', '--numreps',
                        type=int,
                        metavar='INT',
                        default=5,
                        help='''Specify desired number of replicates (sample size)''')
    parser.add_argument('-f', '--outfile',
                        type=str,
                        default='sample',
                        help='''Specify output file name, simulated data will be save to a *.ped file in linkage format''')
    parser.add_argument('--seed',
                        type=float,
                        default=None,
                        help='''Specify seed for random number generator, if left unspecified the current system time will be used''')
    


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
    if not DEBUG:
        pbar = progressbar.ProgressBar(widgets=['Simulating for {} replicates'.format(args.numreps), ' ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA(), ' ', progressbar.FileTransferSpeed()], maxval=int(args.numreps)).start()
    for i in xrange(1, args.numreps+1):
        # per replicate
        
        
        
        samples = []
        diseaseVariantIndices = [weightedRandomIdx(gene1['cumuProbs_dMaf']), weightedRandomIdx(gene2['cumuProbs_dMaf'])]
        for j in xrange(args.samplesize):
            numOffspring = getNumOffspring(offNumProp)
            #### FIXME!
            pedInfo = simPedigree([gene1, gene2], numOffspring, args.mode, args.allelicheteroprop, diseaseVariantIndices)
            samples.append(pedInfo)
        
        
        if not DEBUG:
            pbar.update(i)
    if not DEBUG:
        pbar.finish()    
    #
        # FIXME! here do something about simulated pedigree samples,
        #### FIXME!
        # e.g. save to file, do analysis, etc...
        # saveToFile(...)
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


def simPedigree(genes, numOffspring, mode, hetero, dVarIdx):
    '''
    simulate two-generational pedigree based on given input gene info, number of offspring, mode of inheritance and allelic heterogenity ratio.
    Note: at least two affected offspring are required per family sample
    '''
    causalGeneIdx = 0 if random.random() < hetero[0] else 1
    markerGeneIdx = 1 if causalGeneIdx == 0 else 0
    #
    
    # if mode in ['dominant', 'recessive'], use given dVarIdx as disease variant indices for simulating the first causal haplotype, otherwise re-generate such indices if mode in ['compound_dominant', 'compound_recessive']
    
    
    
    causalHaps = getCausalHaps(genes[causalGeneIdx], mode, numOffspring)
    #markerHaps = getMarkerHaps()
    
    pedInfo = {}
    if mode == 'dominant': # require at least 1 haplotype (out of 4) in parents carries variants, and haps 2,3,4 can carry vars (by maf) only on sites where hap 1 does 
        pass
    elif mode == 'compound_dominant': # at least 1 hap (out of 4) in parents need to carry variants, and haps 2,3,4 can carry vars (by maf) on any causal site
        pass
    elif mode == 'recessive': # at least 1 hap in father and 1 hap in mother have to carry variants on same DSL sites, and haps 3,4 can carry vars only on those sites and at least 3 or 4 has to be wild-type hap
        pass
    elif mode == 'compound_recessive': # at least 1 hap in father and 1 haps in mother have to carry vars on some DSL sites, and haps 3,4 can carry vars on any DSL site and at least 3 or 4 has to be wild-type hap
        pass
    else:
        raise ValueError('')
    
    return pedInfo
  

def getCausalHaps(geneInfo, mode, numOffspring):
    '''
    '''
    causalHaps = collections.OrderedDict({})
    parentalHapCausality = [1,0,0,0] # 1 - causal; 0 - non-causal
    if mode == 'dominant': # require at least 1 haplotype (out of 4) in parents carries variants, and haps 2,3,4 can carry vars (by maf) only on sites where hap 1 does
        while True:
            hap1, dVarIdx = genCausalHap(geneInfo)
            parentalHapCausality[1] = 1 if random.random() < geneInfo['maf'][dVarIdx] else 0
            parentalHapCausality[2] = 1 if (parentalHapCausality.count(1) <  2 and random.random() < geneInfo['maf'][dVarIdx]) else 0
            parentalHapCausality[3] = 1 if (parentalHapCausality.count(1) <  2 and random.random() < geneInfo['maf'][dVarIdx]) else 0
            offspringAff, offHapIdx = genOffspringAffAndHaps(mode, parentalHapCausality, numOffspring)
            if offspringAff.count(2) >= 2:
                # generate hap2, hap3 and hap4
                hap2 = genOtherHap(geneInfo, [dVarIdx]) if parentalHapCausality[1] == 1 else genOtherHap(geneInfo)
                hap3 = genOtherHap(geneInfo, [dVarIdx]) if parentalHapCausality[2] == 1 else genOtherHap(geneInfo)
                hap4 = genOtherHap(geneInfo, [dVarIdx]) if parentalHapCausality[3] == 1 else genOtherHap(geneInfo)
                # fill in parental and offspring causalHaps dict (1 - father; 2 - mother; 3,...,n - offspring)
                causalHaps[1], causalHaps[2] = [hap1, hap2], [hap3, hap4]
                parHaps = [hap1, hap2, hap3, hap4]
                [causalHaps.update({i+3:[parHaps[j[0]], parHaps[j[1]]]}) for i,j in enumerate(offHapIdx)]
                
                break

            
    print dVarIdx
    print parentalHapCausality
    print offspringAff, offHapIdx
    print causalHaps
        
        
    return causalHaps
        

def genOffspringAffAndHaps(mode, parentalHapCausality, numOffspring):
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
            


def genOtherHap(geneInfo, addMutantToSites=[]):
    hap = [0] * len(geneInfo['maf'])
    for idx in geneInfo['nd_idx']:
        if random.random() < geneInfo['maf'][idx]:
            hap[idx] = 1
    # add mutants on given sites
    for idx in addMutantToSites:
        hap[idx] = 1
    return hap
  

def genCausalHap(geneInfo):
    '''
    generate 1st parental haplotype that contains 1 disease variant
    '''
    dVarIdx = weightedRandomIdx(geneInfo['cumuProbs_dMaf'])
    hap = genOtherHap(geneInfo, [dVarIdx])
    return hap, dVarIdx


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
    
    print vars(args)
    
    # calling associated functions
    if DEBUG:
        args.func(args)
    else:
        try:
            args.func(args)
        except Exception as e:
            sys.exit("An ERROR has occured: {}".format(e))

