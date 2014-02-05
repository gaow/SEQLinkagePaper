#!/usr/bin/env python


import argparse, random, tempfile, os, shutil, time, glob, tarfile

DEBUG = True


def arguments(parser):
    parser.add_argument('-l', '--genelength',
                        type=int,
                        metavar='INT',
                        default=1500,
                        help='''Specify length of gene that harbors marker loci, e.g. 1500''')
    parser.add_argument('-o', '--offspring',
                        type=int,
                        metavar='INT, INT',
                        nargs='+',
                        default=[2,5],
                        help='''Specify a range of allowed number of offspring each parent may produce, e.g. 2 5''')
    parser.add_argument('-m', '--markervarfreq',
                        type=float,
                        metavar='FLOAT',
                        default=0.5,
                        help='''Specify variant frequency on marker locus''')
    parser.add_argument('-d', '--diseasevarfreq',
                        type=float,
                        metavar='FLOAT',
                        default=0.5,
                        help='''Specify variant frequency on disease susceptible locus''')
    parser.add_argument('-r', '--recfreq',
                        type=float,
                        metavar='FLOAT',
                        default=0.4,
                        help='''Specify recombination frequency between marker and disease loci (0 <= r <= 0.5)''')
    parser.add_argument('-i', '--mode',
                        type=str,
                        default='recessive',
                        choices=['recessive', 'dominant'],
                        help='''Specify mode of inheritance of disease locus, e.g. recessive or dominant''')
    parser.add_argument('-n', '--numreps',
                        type=int,
                        metavar='INT',
                        default=500,
                        help='''Specify desired number of replicates (sample size)''')
    parser.add_argument('-a', '--ascertain',
                        action='store_true',
                        default=False,
                        help='''Use mlink-calculated LOD scores to ascertain desired family samples''')
    parser.add_argument('--lod',
                        type=float,
                        default=2.0,
                        help='''Specify LOD score threshold used in ascertainment''')
    parser.add_argument('-f', '--outfile',
                        type=str,
                        default='sample',
                        help='''Specify output file name, simulated data will be save to a *.ped file in linkage format''')
    parser.add_argument('-s', '--seed',
                        type=float,
                        default=None,
                        help='''Specify seed for random number generator, if left unspecified the current system time will be used''')



def simLinkage(args):
    '''
    parse cmd input and simulate pedigree samples with two-point linkage
    '''
    # set seed
    if args.seed:
        random.seed(args.seed)
    else:
        random.seed(time.time())
    peds = []
    flag, i = 0, 0
    tempFolder = tempfile.mkdtemp()
    # temporarily change wd to tempFolder
    cwd = os.getcwd()
    os.chdir(tempFolder)
    for i in range(1, args.numreps+1):
        flag = 1
        while flag:
            # create a two-gen pedigree sample and write it to a temp *.ped file
            numOffspring = random.choice(range(args.offspring[0], args.offspring[1]+1))
            ped = createTwoGenPed(mfreq=args.markervarfreq, dfreq=args.diseasevarfreq, numOffspring=numOffspring, recFreq=args.recfreq, modeInheri=args.mode, fam_id=i)
            pedFile = 'ped_'+str(i)+'.ped'
            writePedsToFile(peds=[ped], fileName=pedFile)
            # compute lod score and ascertain, if successful save to peds
            ##!!FIXME!!
            if args.ascertain:
                if calLOD(pedFile, threshold=args.lod):
                    flag = 0
            else:
                flag = 0
            if not flag:
                peds.append(ped)
    ##
    # write peds to output file (change wd back to cwd here)
    os.chdir(cwd)
    writePedsToFile(peds=peds, fileName=args.outfile)
    # compress ped files of single families to args.outfile.bz2, and delete tempFolder
    bz2Save(args.outfile, tempFolder)
    return 


def calLOD(pedFile, threshold):
    '''
    calculate LOD score of a given pedigree and return True if
    ascertaining criterion is satisfied
    '''
    # !!FIXME!! #
    return True


def writePedsToFile(peds, fileName):
    '''
    write pedigree samples to 'fileName' in linkage format
    '''
    if not fileName.endswith('.ped'):
        fileName += '.ped'
    fi = open(fileName, 'w')
    for ped in peds:
        for ind_info in ped:
            fi.write(' '.join(str(x) for x in ind_info) + '\n')
    fi.close()
    return


def createTwoGenPed(mfreq, dfreq, numOffspring, recFreq, modeInheri, fam_id):
    '''
    create a two-gen pedigree and return a list of inds' info in the order
    of [father_info, mother_info, offspring_1_info, offspring_2_info,...,offspring_n_info],
    where each ind_info element corresponds to info of each row in linkage(ped) file
    '''
    ped = []
    # father & mother info
    fa_haps, ma_haps = getFounderHaps(mfreq, dfreq), getFounderHaps(mfreq, dfreq)
    fa_pheno, ma_pheno = getAffectionStatus(haps=fa_haps, mode=modeInheri), getAffectionStatus(haps=ma_haps, mode=modeInheri)
    fa_info, ma_info = [fam_id, 1, 0, 0, 1, fa_pheno], [fam_id, 2, 0, 0, 2, ma_pheno]
    fa_info.extend(fa_haps), ma_info.extend(ma_haps)
    ped.append(fa_info), ped.append(ma_info)
    # offspring info
    for ind_id in range(3, numOffspring+3):
        ind_haps = getOffspringHaps(fa_haps, ma_haps, recFreq)
        ind_pheno = getAffectionStatus(haps=ind_haps, mode=modeInheri)
        ind_sex = 1 if random.random() <= 0.5 else 2
        ind_info = [fam_id, ind_id, 1, 2, ind_sex, ind_pheno]
        ind_info.extend(ind_haps)
        ped.append(ind_info)
    return ped


def getOffspringHaps(fa_haps, ma_haps, recFreq):
    '''
    '''
    def _getGamete(haps, recFreq):
        # if recombine
        if random.random() < recFreq:
            gamete1 = [haps[0], haps[3]]
            gamete2 = [haps[1], haps[2]]
        else:
            gamete1 = [haps[0], haps[2]]
            gamete2 = [haps[1], haps[3]]
        return gamete1 if random.random() < 0.5 else gamete2
    #
    gam1, gam2 = _getGamete(fa_haps, recFreq), _getGamete(ma_haps, recFreq)
    return [gam1[0], gam2[0], gam1[1], gam2[1]]


def getFounderHaps(mfreq, dfreq):
    '''
    mfreq - variant frequency on marker locus
    dfreq - variant frequency on disease locus
    '''
    return [1 if random.random() < mfreq else 0,
            1 if random.random() < mfreq else 0,
            1 if random.random() < dfreq else 0,
            1 if random.random() < dfreq else 0]


def getAffectionStatus(haps, mode):
    '''
    currently assume 100% penetrance for two-point linkage analysis
    mode: dominant, recessive
    affection status: 1-unaffected, 2-affected
    '''
    dsl = haps[-2:]
    if mode == 'recessive' and dsl.count(1) == 2:
        return 2
    elif mode == 'dominant' and dsl.count(1) > 0:
        return 2
    else:
        return 1
    
    
def bz2Save(fileName, tempFolder):
    '''
    zip 'rep#' files in tempFolder to fileName and delete tempFolder
    '''
    repNames = glob.glob(os.path.join(tempFolder, '*'))
    tar = tarfile.open(fileName+'.bz2', 'w:bz2')
    cwd = os.getcwd()
    os.chdir(tempFolder)
    for name in repNames:
        #tar.add(name)
        tar.add(os.path.basename(name))
    # remove tempFolder
    shutil.rmtree(tempFolder)
    os.chdir(cwd)
    return



if __name__ == '__main__':
    master_parser = argparse.ArgumentParser(
        description = '''Program to generate two generational family samples with two loci under linkage and ascertain for families that may achieve desired LOD scores''',
        prog = 'simlinkage',
        epilog = '''Biao Li (biaol@bcm.edu) (c) 2014.'''
    )
    master_parser.add_argument('--version,', action='version', version='%(prog)s 0.1.0')
    arguments(master_parser)
    master_parser.set_defaults(func=simLinkage)
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






