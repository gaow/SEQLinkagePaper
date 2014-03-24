#!/usr/bin/python
# Di Zhang diz@bcm.edu
#

from SEQLinco.Utils import *
from collections import deque, defaultdict
from multiprocessing import Queue, Process, cpu_count
from os.path import splitext, basename, isdir
from itertools import chain
from shutil import copyfile, rmtree
from scipy.optimize import minimize_scalar
#utils that allow lambda function in mutilprocessing map
#grabbed from http://stackoverflow.com/a/16071616 by klaus-se
def spawn(f):
    def fun(q_in,q_out):
        while True:
            i,x = q_in.get()
            if i is None:
                break
            q_out.put((i,f(x)))
    return fun
def parmap(f, X, nprocs = cpu_count()):
    q_in   = Queue(1)
    q_out  = Queue()
    proc = [Process(target=spawn(f),args=(q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]
    [p.join() for p in proc]
    return [x for i,x in sorted(res)]

#formatters
#the handler, called from main, can call specific formmater.
def format_linkage(tpeds, tfam, prev, wild_pen, muta_pen,out_format='mlink', inherit_mode='AR'):
    #pool = Pool(env.jobs)
    if out_format == 'plink':
        parmap(lambda x: format_plink(x, tfam), tpeds, env.jobs)
    elif out_format == 'mlink':
        parmap(lambda x: format_mlink(x, tfam, prev, wild_pen, muta_pen, inherit_mode), tpeds, env.jobs)

#plink format, ped and map 
def format_plink(tped, tfam):
    out_base = 'PLINK/{}'.format(splitext(basename(tped))[0])
    with open(tped) as tped_fh, open(tfam) as tfam_fh:
        geno = []
        with open(out_base + '.map', 'w') as m:
            for line in tped_fh:
                s = line.strip().split()
                m.write(env.delimiter.join(s[:3]))
                geno.append(deque(s[4:]))
        m.close()
        with open(out_base + '.ped', 'w') as p:
            for line in tfam_fh:
                p.write(line.strip())
                map(lambda x: p.write(' {} {}'.format(x.popleft(), x.popleft)), geno)
                p.write("\n")
        p.close()
    tped_fh.close()
    tfam_fh.close()

#mlink format, .pre and .loc
#per locus, per family based
#because the haplotype patterns are different from family to family.
#You can analyze them all together
def format_mlink(tped, tfam, prev, wild_pen, muta_pen, inherit_mode):
    out_base = 'MLINK/{}'.format(env.tmp_dir, splitext(basename(tped))[0])
    try:
        rmtree(out_base)
    except:
        pass
    env.log("Start converting to mlink format for {} ...\n".format(basename(out_base)), flush=True)
    with open(tped) as tped_fh, open(tfam) as tfam_fh:
        fams = parse_tfam(tfam_fh)
        #parse per family per locus AF file
        af = defaultdict(lambda: [])
        #try to open the file for allele frequencies, otherwise use the defaut value
        try:
            with open('cache/{}.freq'.format(basename(out_base))) as af_fh:
                for line in af_fh:
                    s = line.strip().split()
                    af[(s[0],s[1])] = s[2:]
        except IOError:
            env.error('freq info not properly read for [{}]'.format(basename(out_base)))
        #parse tped
        heter_pen = wild_pen
        if inherit_mode == 'AD':
            heter_pen = muta_pen
        #else:
        #    env.error('Inheritance mode {} not implemented yet'.format(inherit_mode))
        for line in tped_fh:
            s = line.strip().split()
            gene, gno = re.search(r'^(\S+?)(?:\[(\d+)\])?$', s[1]).groups()
            if not gno:
                gno = '0'
            for fid in fams:
                workdir = '{}/{}/{}'.format(out_base, gene, fid)
                mkpath(workdir)
                #env.error("fid {} num {}\n".format(fid, fams[fid].get_member_ids()))
                fam_af = af[(fid, s[1])]
                if not fam_af:
                    continue
                    #env.error('no family-wise allele freq info for {} {} {}'.format(fid, s[1], fam_af))
                    #fam_af = ['0.1'] * gs_num
                ids = fams[fid].get_sorted_ids()
                idxes = map(lambda x: fams[fid].get_member_idx(x), ids)
                gs = map(lambda x: s[2 * x + 4 : 2 * x + 6], idxes)
                #if re.search(r'\[21\]$', workdir):
                #    env.log(set(filter(lambda x: x != '0', chain(*gs))))
                gs_num = len(set(filter(lambda x: x != '0', chain(*gs))))
                if gs_num >= 10:
                    env.log('pattern number larger than 9 for unit {} in family {}, skipped'.format(s[1], fid))
                    continue
                with open('{}/{}.PRE'.format(workdir, gno), 'w') as pre:
                    pre.write(''.join("{} {} {} {}\n".format(fid, fams[fid].print_member(pid), s[2*fams[fid].get_member_idx(pid) + 4], s[2*fams[fid].get_member_idx(pid) + 5]) for pid in ids))
                with open('{}/{}.LOC'.format(workdir, gno), 'w') as loc:
                    loc.write("2 0 0 5\n")
                    loc.write("0 0.0 0.0 0\n")
                    loc.write("1 2\n")
                    loc.write("1 2\n")
                    loc.write(" {} {}\n".format(1 - prev, prev))
                    loc.write(" 1\n")
                    loc.write(" {} {} {}\n".format(wild_pen, heter_pen, muta_pen))
                    loc.write("3 {}\n".format(gs_num))
                    loc.write(' ' + ' '.join(fam_af) + "\n")
                    loc.write("0 0\n")
                    loc.write("0.0\n")
                    loc.write("1 0.05 0.45\n")
                if not os.listdir(workdir):
                    env.log('empty dir [{}], deteled'.format(workdir))
                    os.rmdir(workdir)
            if not os.listdir('{}/{}'.format(out_base, gene)):
                env.log('empty dir [{}/{}], deteled'.format(out_base, gene))
                os.rmdir(workdir)
    tped_fh.close()
    tfam_fh.close()
    env.log("Finished mlink format for {}.\n".format(out_base), flush=True)
  
#parse tfam file, store families into the Pedigree class                
def parse_tfam(fh):
    fams = defaultdict(lambda: Pedigree())
    idx = 0
    for line in fh:
        s = line.strip().split()
        fams[s[0]].add_member(s[1:], idx)
        idx += 1
    return fams

#This is to sort the members of a pedigree.
#To make sure that parents come before offsprings.
class Pedigree:
    def __init__(self):
        self.fid = None
        self.data = {}
        self.graph = defaultdict(lambda:[]) #, defaultdict(lambda:[])]
        self.sorted = []
        
    def add_member(self, info, idx): #list [pid, father, mother, sex, pheno]
        if info[1] != '0' and info[2] != '0':
            self.graph[info[1]].append(info[0])
            self.graph[info[2]].append(info[0])
        self.data[info[0]] = info + [idx]
        
    def get_member_info(self, pid):
        return self.data[pid][:-1]

    def get_member_idx(self, pid):
        return self.data[pid][-1]
    
    def get_member_ids(self):
        return self.data.keys()

    def print_member(self, pid):
        return ' '.join(self.get_member_info(pid))

    def get_sorted_ids(self):
        if self.sorted:
            return self.sorted
        else:
            #This algorithm was first described by Kahn (1962)
            S_no_parents = filter(lambda x: True if self.get_member_info(x)[1] == '0' else False, self.get_member_ids())
            graph = self.graph.copy()
            while(S_no_parents):
                n = S_no_parents.pop()
                self.sorted.append(n)
                if n not in graph:
                    continue
                offsprings = graph.pop(n)
                for m in offsprings:
                    father = self.get_member_info(m)[1]
                    mother = self.get_member_info(m)[2]
                    if father not in graph and mother not in graph:
                        S_no_parents.append(m)
            if graph:
                raise Exception("There is a loop in the pedigree: {}\n".format(' '.join(graph.keys())))
            else:
                return self.sorted

#runners
#context manager
#grabbed from http://stackoverflow.com/a/13197763, by Brian M. Hunt
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def run_linkage(runner, blueprint):
    #if 'mlink' not in env.formats:
    #    formatMlink()
    #chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    #cmds = ['runMlink.pl MLINK/{}.{} {} {}'.format(env.output, chrs[i], env.resource_dir, blueprint) for i in range(25)]
    #runCommands(cmds, max(min(env.jobs, cmds), 1))
    if runner == 'mlink':
        workdirs = glob.glob('MLINK/{}.chr*'.format(env.output))
        parmap(lambda x: run_mlink(x, blueprint) , workdirs, env.jobs)
    
def run_mlink(workdir, blueprint):
    env.log("Start running mlink for {} ...".format(workdir), flush=True)
    #hash genes into genemap
    genemap = {}
    with open(blueprint) as f:
        for line in f.readlines():
            chrID, start, end, gene = line.strip().split()[:4]
            genemap[gene] = [chrID, int(start), int(end)]
    #This is a hack, should be simpler. Because some genes might have their suffix [\d+],
    #We couldn't find GENENAME[\d+] in the genemap, we have to look up the original name,
    #and when we sort them, suffixes if exist should be added properly.
    PNULL = open(os.devnull, 'w')
    mkpath('MLINK/heatmap')
    lods_fh = open('MLINK/heatmap/{}.lods'.format(basename(workdir)), 'w')
    hlods_fh = open('MLINK/heatmap/{}.hlods'.format(basename(workdir)), 'w')
    #for unitdir in sorted(filter(isdir, glob.glob(workdir + '/*')), key=lambda x: genemap[re.sub(r'^(\S+?)(?:\[\d+\])?$', r'\1', basename(x))] + [int(re.sub(r'^(?:\S+?)(?:\[(\d+)\])?$', r'\1', basename(x)) if re.search(r'\]$', x) else 0)]):
    genes = filter(lambda g: g in genemap, map(basename, glob.glob(workdir + '/*')))
    for gene in sorted(genes, key=lambda g: genemap[g]):
        lods = {}
        hlods = {}
        fams = map(basename, filter(isdir, glob.glob('{}/{}/*'.format(workdir, gene))))
        for fam in fams:
            with cd('{}/{}/{}'.format(workdir, gene, fam)):
                units = map(lambda x: re.sub(r'^(\d+?)\.PRE$', r'\1', x) ,glob.glob('*.PRE'))
                for unit in units:
                    copyfile('{}.LOC'.format(unit), 'datafile.dat')
                    copyfile('{}.PRE'.format(unit), 'pedfile.pre')
                    subprocess.call(['makeped', 'pedfile.pre', 'pedfile.ped', 'n'], stdout=PNULL, stderr=PNULL)
                    subprocess.call(['pedcheck', '-p', 'pedfile.ped', '-d', 'datafile.dat', '-c'], stdout=PNULL, stderr=PNULL)
                    copyfile('zeroout.dat', 'pedfile.dat')
                    runCommand('unknown')
                    runCommand('mlink')
                    copyfile('outfile.dat', '{}.out'.format(unit))        
                    #clean mlink tmp files
                    for f in set(glob.glob('*.dat') + glob.glob('ped*') + ['names.tmp']):
                        os.remove(f)
                    #collect lod scores of different thelta for the fam
                    with open('{}.out'.format(unit)) as out:
                        raw = out.read()
                        for i in re.finditer(r'^THETAS\s+(0\.\d+)(?:\n.+?){7}LOD SCORE =\s+(-?\d+\.\d+)', raw, re.MULTILINE):
                            theta, lod = map(float, i.group(1,2))
                            if theta not in lods:
                                lods[theta] = {fam: lod}
                            elif fam not in lods[theta] or lod > lods[theta][fam]:
                                lods[theta][fam] = lod
        for theta in sorted(lods.keys()):
            res = minimize_scalar(hlod_fun(lods[theta].values(), -1), bounds=(0,1), method='bounded', options={'xtol':1e-8})
            a = res.x
            lods_fh.write('{} {} {}\n'.format(gene, theta, sum(lods[theta].values())))
            hlods_fh.write('{} {} {} {}\n'.format(gene, theta, a, hlod_fun(lods[theta].values())(a)))
    PNULL.close()
    lods_fh.close()
    hlods_fh.close()
    heatmap('MLINK/heatmap/{}.lods'.format(basename(workdir)))
    heatmap('MLINK/heatmap/{}.hlods'.format(basename(workdir)))
    env.log("Finished running mlink for {}.".format(workdir), flush=True)
    
def heatmap(file):
    env.log("Start ploting heatmap for {} ...\n".format(file), flush=True)
    lods = []
    with open(file, 'r') as f:
        for line in f.readlines():
            lod = line.split()[-1]
            lods.append(lod)
        lods = np.array(map(float,lods)).reshape((11,-1))
        ppl.pcolormesh(lods)
        plt.savefig('{}.png'.format(file))
        plt.close()
    env.log("Finished ploting heatmap for {}.\n".format(file), flush=True)

def plotMlink():
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    dirs = filter(lambda x: os.path.exists(x), ['MLINK/{}.{}'.format(env.output, i) for i in chrs])
    mkpath('MLINK/heatmap')
    pool = Pool(env.jobs)
    pool.map(heatmap, dirs)

def hlod_fun(Li, sign=1):
    def _fun(alpha):
        return sign * sum(np.log10(alpha*np.power(10, Li) + 1 - alpha))
    return _fun 
