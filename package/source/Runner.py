#!/usr/bin/python
# Di Zhang diz@bcm.edu
#

from SEQLinco.Utils import *
from collections import deque, defaultdict
from multiprocessing import Queue, Process, cpu_count
from os.path import splitext, basename
from itertools import chain
#utils that allow lambda function in mutilprocessing map
#copied from http://stackoverflow.com/a/16071616 by klaus-se on StackOverflow.com
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
def format_linkage(tpeds, tfam, out_format='mlink', prev=0.001, inherit_mode='AR'):
    #pool = Pool(env.jobs)
    if out_format == 'plink':
        parmap(lambda x: format_plink(x, tfam), tpeds, env.jobs)
    elif out_format == 'mlink':
        parmap(lambda x: format_mlink(x, tfam, prev, inherit_mode), tpeds, env.jobs)
               
def format_plink(tped, tfam):
    out_base = 'PLINK/' + splitext(basename(tped))[0]
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

def format_mlink(tped, tfam, prev, inherit_mode):
    out_base = 'MLINK/' + splitext(basename(tped))[0]
    with open(tped) as tped_fh, open(tfam) as tfam_fh:
        fams = parse_tfam(tfam_fh)
        #parse per family per locus AF file
        af = defaultdict(lambda: [])
        #try to open the file for allele frequencies, otherwise use the defaut value
        try:
            with open('cache/variant_af.tbl') as af_fh:
                for line in af_fh:
                    s = line.strip().split()
                    af[(s[0],s[1])] = s[2:]
        except IOError:
            pass
        #parse tped
        heter_pen = 0.0
        if inherit_mode == 'AD':
            heter_pen = 0.9
        #else:
        #    env.error('Inheritance mode {} not implemented yet'.format(inherit_mode))
        for line in tped_fh:
            s = line.strip().split()
            workdir = '{}/{}'.format(out_base, s[1])
            mkpath(workdir)
            for fid in fams:
                #env.error("fid {} num {}\n".format(fid, fams[fid].get_member_ids()))
                fam_af = af[(fid, s[1])]
                gs_num = 0
                with open('{}/{}.PRE'.format(workdir, fid), 'w') as pre:
                    ids = fams[fid].get_sorted_ids()
                    idxes = map(lambda x: fams[fid].get_member_idx(x), ids)
                    gs = map(lambda x: s[2 * x + 4 : 2 * x + 6], idxes)
                    gs_num = len(set(filter(lambda x: x != 0, chain(*gs))))
                    if not fam_af:
                        fam_af = ['0.1'] * gs_num
                    pre.write(''.join("{} {} {} {}\n".format(fid, fams[fid].print_member(pid), s[2*fams[fid].get_member_idx(pid) + 4], s[2*fams[fid].get_member_idx(pid) + 5]) for pid in ids))
                with open('{}/{}.LOC'.format(workdir, fid), 'w') as loc:
                    loc.write("2 0 0 5\n")
                    loc.write("1 2\n")
                    loc.write("1 2\n")
                    loc.write(" {} {}\n".format(1 - prev, prev))
                    loc.write(" 1\n")
                    loc.write(" 0.0 {} 0.9\n".format(heter_pen))
                    loc.write("3 {}\n".format(gs_num))
                    loc.write(' ' + ' '.join(fam_af) + "\n")
                    loc.write("0 0\n")
                    loc.write("0.0\n")
                    loc.write("1 0.05 0.45\n")
    tped_fh.close()
    tfam_fh.close()


    
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
def runMlink(blueprint):
    #if 'mlink' not in env.formats:
    #    formatMlink()
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    cmds = ['runMlink.pl MLINK/{}.{} {} {}'.format(env.output, chrs[i], env.resource_dir, blueprint) for i in range(25)]
    runCommands(cmds, max(min(env.jobs, cmds), 1))

def heatmap(dir):
    env.log("start ploting heatmap for" + dir + "\n")
    lods = []
    with open(dir + '/all_lodscores.txt', 'r') as f:
        for line in f.readlines():
            lod = line.split()[-1]
            lods.append(lod)
        lods = np.array(map(float,lods)).reshape((10,-1))
        ppl.pcolormesh(lods)
        plt.savefig('MLINK/heatmap/{}.lods.png'.format(os.path.basename(dir)))
    env.log("end ploting heatmap for" + dir + "\n")

def plotMlink():
    chrs = ['chr{}'.format(i+1) for i in range(22)] + ['chrX', 'chrY', 'chrXY']
    dirs = filter(lambda x: os.path.exists(x), ['MLINK/{}.{}'.format(env.output, i) for i in chrs])
    mkpath('MLINK/heatmap')
    pool = Pool(env.jobs)
    pool.map(heatmap, dirs)
