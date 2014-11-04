#!/usr/bin/python

import sys, os, argparse, gzip, subprocess

##Global conf of each running instance.
##These are just some naming conventions of files and the genes to be
##used. Could be changed when another analysis being performed.
##Other parts of this script should be reusable.
class Conf:
    def __init__(self, args):
        self.args = args
        ##path and program conventions
        self.srcDir = '../0data/rep{}'.format(args.rep)
        self.tgtDir = 'rep{}'.format(args.rep)
        self.Program = 'gh'
        self.srcFile = 'FAM{}RUN{}rep{}.{{}}'.format(args.fam, args.run, args.rep)
        self.tgtFile = 'FAM{}RUN{}rep{}gene{{}}.{{}}'.format(args.fam, args.run, args.rep)
        #VCF parameters
        self.genes = [0, 1] #There are two genes
        self.Chr = [7, 13] #the genes to parse are on chr7 and chr13
        ##Ped parameters
        self.PedOffset = 6 #the first 6 columns of a ped file are fam info
        ##Loci parameters
        self.moi = 'AR'
        self.moiPrev = {'AR': '0.001 0.001 0.999',
                        'AD': '0.001 0.999 0.999'}
        self.markAF = 0.01
        self.markRR = 0.00001
        self.locTpl = self.getLocTemplate()
        ##Genehunter parameters
        self.step = 1
        self.mapFunc = 'haldane'
        
    def getLocTemplate(self):
        head = '{} 0 0 5\n0 0.0 0.0 0\n'
        loci = '{}\n1 2\n0.99 0.01\n'
        disloc = '1\n{}\n'.format(self.moiPrev[self.moi])
        markers = '{}'
        sex = '0 0\n'
        recRate = '0.0 {}\n'
        recInc = '1 0.05 0.45\n'
        Tpl = head + loci + disloc + markers + sex + recRate + recInc
        return Tpl


##Utils
#Unfortunately, python lacks the arbitary array slice
def getSlice(lst, idxes):
    res = []
    #print idxes
    for i in idxes:
        res.append(lst[i])
    return res

#get source file name by ext
def getSrc(cnf, ext):
    src = os.path.join(cnf.srcDir, cnf.srcFile.format(ext))
    return src

#get target file name by ext and gene No.
def getTgt(cnf, ext, gene):
    tgt = os.path.join(cnf.tgtDir, cnf.tgtFile.format(gene, ext))
    return tgt

#cd
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

##Classes
class Loc:
    def __init__(self, cnf, num, gene):
        self.moiPrev = cnf.moiPrev[cnf.moi] 
        self.markAF = cnf.markAF
        self.markRR = cnf.markRR
        self.tpl = cnf.locTpl
        self.File = getTgt(cnf, 'loc', gene)
        self.numLoci = num + 1#number of markers
        self.numMarkers = num
       
    def getLoci(self):
        return ' '.join(map(lambda x: str(1 + x), range(self.numLoci)))
    
    def getMarkers(self):
        res = ''
        for m in range(self.numMarkers):
            res += '3 2\n{} {}\n'.format(1 - self.markAF, self.markAF)
        return res
        
    def getRates(self):
        return ' '.join([str(self.markRR)] * (self.numMarkers - 1))

    def make(self):
        #print self.File
        with open(self.File, 'w') as l:
            l.write(self.tpl.format(self.numLoci, self.getLoci(), self.getMarkers(), self.getRates()))
 
class Pre:
    def __init__(self, cnf, idxes, gene):
        self.ped = getSrc(cnf, 'ped')
        self.File = getTgt(cnf, 'pre', gene)
        self.offset = cnf.PedOffset
        self.idxes = idxes

    def make(self):
        with open(self.ped) as p, open(self.File, 'w') as f:
            for line in p:
                s = line.strip('\n').split()
                fam = map(lambda x: x.split(':')[-1], s[:self.offset])
                geno = map(lambda x: str(1 + int(x)), getSlice(s[self.offset:], self.idxes))
                f.write('{} {}\n'.format(' '.join(fam), ' '.join(geno)))

##Meaningful funtions
#parse the vcf to determine the non-all-wildtype sites of a gene
def FindSites(cnf, gene):
    pos = cnf.Chr[gene]
    vcf = getSrc(cnf, 'vcf.gz')
    res = []
    idx = 0
    with gzip.open(vcf) as v:
        for l in v:
            if l.startswith('#'):
                continue
            s = l.strip('\n').split()
            if s[0] != str(pos):
                idx += 1
                continue
            g = s[9:]
            if len(filter(lambda x: x != '0/0', g)) > 0:
                res.append(idx * 2)
                res.append(idx * 2 + 1)
            idx += 1
    #print res
    return res

#prepare file for genehunter, most of the job was done by the Loc and
#Pre classes.
def PrepFiles(cnf, idxes, gene):
    l = Loc(cnf, len(idxes)/2, gene)
    l.make()
    p = Pre(cnf, idxes, gene)
    p.make()
    return (l, p)

##Find the highest lod score and store it in a file
def FindLod(res, lods):
    start, end, val = (0, 0, -100)
    with open(lods, 'w') as o:
        for idx, line in enumerate(res.split('\n')):
            #print idx
            #print '{} {} {}'.format(start, end, line[:6])
            if line[:6] == 'npl:7>':
                o.write('{}\n'.format(val))
                break
            if start == 1:
                if line.strip() != '':
                    lod = line.strip().split()[1]                    
                    if lod > val:
                        val = lod
            if line[:6] == 'npl:6>':
                start = 1
                #print start
 
#Actually run the genehunter
def RunGeneHunter(Loc, Pre, cnf):
    with cd(cnf.tgtDir):
        cmds = 'load {}\nincrement steps {}\nmap function {}\nanalysis lod\nscan {}\ntotal stat\n'
        cmds = cmds.format(os.path.basename(Loc.File), cnf.step, cnf.mapFunc, os.path.basename(Pre.File))
        #print cmds
        gh = subprocess.Popen('gh', stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        res = gh.communicate(cmds)[0]
        lods = os.path.splitext(os.path.basename(Loc.File))[0] + '.lods'
        FindLod(res, lods)
        try:
            gh.kill()
        except:
            pass

def main():
    parser = argparse.ArgumentParser(description='fam run rep')
    parser.add_argument('--fam',
                        help='fam No.')
    parser.add_argument('--run',
                        help='run No.')
    parser.add_argument('--rep',
                        help='rep No.')
    args = parser.parse_args()
    cnf = Conf(args)
    for gene in cnf.genes:
        Idxes = FindSites(cnf, gene)
        Loc, Pre = PrepFiles(cnf, Idxes, gene)
        RunGeneHunter(Loc, Pre, cnf)


if __name__ == '__main__':
    main()
