#!/usr/bin/env python


import argparse, random, tempfile, os, shutil, time, glob, tarfile

DEBUG = True

OFF_PROP = {'2': 0.6314, '3': 0.2556, '4_more': 0.1130}

def arguments(parser):
    parser.add_argument('-g', '--genes',
                        type=str,
                        metavar='FILE, FILE',
                        nargs='+',
                        default=[],
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
                        choices=['recessive', 'dominant', 'compound_recessive'],
                        help='''Specify mode of inheritance of disease locus, e.g. recessive or dominant or compound_recessive or compound_dominant''')
    parser.add_argument('-s', '--samplesize',
                        type=int,
                        metavar='INT',
                        default=3000,
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
                        default=500,
                        help='''Specify desired number of replicates (sample size)''')
    parser.add_argument('-f', '--outfile',
                        type=str,
                        default='sample',
                        help='''Specify output file name, simulated data will be save to a *.ped file in linkage format''')
    parser.add_argument('-s', '--seed',
                        type=float,
                        default=None,
                        help='''Specify seed for random number generator, if left unspecified the current system time will be used''')



def simSEQLinco(args):
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
    for i in range(1, args.numreps+1):
        # per replicate
        samples = []
        for j in range(args.samplesize):
            numOffspring = getNumOffspring(offNumProp)
            pedInfo = simPedigree([gene1, gene2], numOffspring, args.mode, args.allelicheteroprop)
            samples.append(pedInfo)
        #
        # FIXME! here do something about simulated pedigree samples,
        # e.g. save to file, do analysis, etc...
    
    return
    
    
def checkInput(args):
    '''
    validate user inputs and raise InputError in case of exception
    '''
    return


def updateOffNumProp(offspringRange):
    '''
    return an updated proportion of number of offspring, according to both OFF_PROP and args.offspring
    '''

    return offProp
    


def parseGeneInfo(fileName):
    '''
    parse input gene file (*.tsv) and return dict obj of gene info
    '''
    info = {'pos':[], 'maf':[], 'function':[]}
    ## FIXME!
    return info


def simPedigree(genes, numOffspring, mode, hetero):
    '''
    simulate two-generational pedigree based on given input gene info, number of offspring, mode of inheritance and allelic heterogenity ratio.
    Note: at least two affected offspring are required per family sample
    '''
    varGeneIdx = 0 if random.random() < hetero[0] else 1
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
    

def getNumOffspring(offNumProp):
    '''
    return randomly generated number of offspring according to proportion of number of offspring
    '''
    return num



if __name__ == '__main__':
    master_parser = argparse.ArgumentParser(
        description = '''Program to generate two generational family samples with two gene regions''',
        prog = 'simSEQLinco',
        epilog = '''Biao Li (biaol@bcm.edu) and Hang Dai (hang.dai@bcm.edu) (c) 2014.'''
    )
    master_parser.add_argument('--version,', action='version', version='%(prog)s 0.1.0')
    arguments(master_parser)
    master_parser.set_defaults(func=simSEQLinco)
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

