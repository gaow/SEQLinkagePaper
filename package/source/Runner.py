#!/usr/bin/python
# Di Zhang diz@bcm.edu
#
from SEQLinkage.Utils import *
from collections import deque, defaultdict
from multiprocessing import Queue, Process, cpu_count
from os.path import splitext, basename, isdir, isfile
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
def format(tpeds, tfam, prev, wild_pen, muta_pen, out_format, inherit_mode, theta_max, theta_inc):
    if out_format == 'plink':
        parmap(lambda x: format_plink(x, tfam), tpeds, env.jobs)
    elif out_format == 'linkage':
        parmap(lambda x: format_linkage(x, tfam, prev, wild_pen, muta_pen, inherit_mode, theta_max, theta_inc), tpeds, env.jobs)

#plink format, ped and map 
def format_plink(tped, tfam):
    out_base = '{}/PLINK/{}'.format(env.tmp_dir, splitext(basename(tped))[0])
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

#linkage format, .pre and .loc
#per locus, per family based
#because the haplotype patterns are different from family to family.
#You can analyze them all together
def format_linkage(tped, tfam, prev, wild_pen, muta_pen, inherit_mode, theta_max, theta_inc):
    out_base = '{}/LINKAGE/{}'.format(env.tmp_dir, splitext(basename(tped))[0])
    try:
        rmtree(out_base)
    except:
        pass
    env.log("Start converting to LINKAGE format for {} ...\n".format(basename(out_base)), flush=True)
    with open(tped) as tped_fh, open(tfam) as tfam_fh:
        fams = parse_tfam(tfam_fh)
        #parse per family per locus AF file
        af = defaultdict(lambda: [])
        #try to open the file for allele frequencies, otherwise use the defaut value
        try:
            with open(os.path.join(env.tmp_dir, 'CACHE', basename(out_base) + '.freq')) as af_fh:
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
                    loc.write("1 {} {}\n".format(theta_inc, theta_max))
                if not os.listdir(workdir):
                    env.log('empty dir [{}], deteled'.format(workdir))
                    os.rmdir(workdir)
            if not os.listdir('{}/{}'.format(out_base, gene)):
                env.log('empty dir [{}/{}], deteled'.format(out_base, gene))
                os.rmdir(workdir)
    tped_fh.close()
    tfam_fh.close()
    env.log("Finished LINKAGE format for {}.\n".format(out_base), flush=True)
  
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

def run_linkage(blueprint, theta_inc, theta_max):
    workdirs = glob.glob('{}/LINKAGE/{}.chr*'.format(env.tmp_dir, env.output))
    parmap(lambda x: linkage_worker(blueprint, x, theta_inc, theta_max) , workdirs, env.jobs)
    
def linkage_worker(blueprint, workdir, theta_inc, theta_max):
    env.log("Start running LINKAGE for {} ...".format(workdir), flush=True)
    #hash genes into genemap
    genemap = {}
    if blueprint:
        with open(blueprint) as f:
            for line in f.readlines():
                chrID, start, end, gene = line.strip().split()[:4]
            genemap[gene] = [chrID, int(start), int(end)]
    else:
        tped = os.path.join(env.tmp_dir, basename(workdir))
        with open(tped) as f:
            for line in f.readlines():
                items = line.strip().split()[:4]
                chrID = items[0]
                gene = items[1]
                pos = items[3]
                genemap[gene] = [chrID, int(pos), int(pos+1)]
    mkpath('{}/heatmap'.format(env.output))
    lods_fh = open('{}/heatmap/{}.lods'.format(env.output, basename(workdir)), 'w')
    hlods_fh = open('{}/heatmap/{}.hlods'.format(env.output, basename(workdir)), 'w')
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
                    runCommand(['makeped', 'pedfile.pre', 'pedfile.ped', 'n'],
                               show_stderr = False, return_zero = False)
                    runCommand(['pedcheck', '-p', 'pedfile.ped', '-d', 'datafile.dat', '-c'],
                               show_stderr = False, return_zero = False)
                    copyfile('zeroout.dat', 'pedfile.dat')
                    runCommand('unknown')
                    runCommand('mlink')
                    copyfile('outfile.dat', '{}.out'.format(unit))        
                    #clean linkage tmp files
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
            lods_fh.write('{} {} {} {}\n'.format(gene, ' '.join(map(str, genemap[gene])), theta, sum(lods[theta].values())))
            hlods_fh.write('{} {} {} {} {}\n'.format(gene, ' '.join(map(str, genemap[gene])), a, theta, hlod_fun(lods[theta].values())(a)))
    lods_fh.close()
    hlods_fh.close()
    heatmap('{}/heatmap/{}.lods'.format(env.output, basename(workdir)), theta_inc, theta_max)
    heatmap('{}/heatmap/{}.hlods'.format(env.output, basename(workdir)), theta_inc, theta_max)
    env.log("Finished running LINKAGE for {}.".format(workdir), flush=True)
    
def heatmap(file, theta_inc, theta_max):
    env.log("Start ploting heatmap for {} ...\n".format(file), flush=True)
    lods = []
    with open(file, 'r') as f:
        for line in f.readlines():
            theta,lod = line.split()[-2:]
            if float(theta) >= theta_max:
                continue
            lods.append(lod)
        Num=int(round(theta_max/theta_inc))
        lods = np.array(map(float,lods)).reshape((-1,Num))
        chrID = re.search(r'\.chr([0-9XY]+)\.', file).group(1)
        fig, ax = plt.subplots(1)
        ax.set_title('Chromosome {}'.format(chrID))
        ppl.pcolormesh(fig,ax,lods.transpose(),
                       yticklabels=np.round(np.array(range(Num)) * theta_inc,2).tolist())
        fig.savefig('{}.png'.format(file))
        #fig.close()
    env.log("Finished ploting heatmap for {}.\n".format(file), flush=True)

def hlod_fun(Li, sign=1):
    def _fun(alpha):
        return sign * sum(np.log10(alpha*np.power(10, Li) + 1 - alpha))
    return _fun 
    
def html(theta_inc, theta_max, limit):
    head = """<html>
    <head>
    <title>Results for {}</title>""".format(env.output) + """
    <style>
    table,th,td
    {
    border:1px solid black;
    border-collapse:collapse;
    }
    </style>
    <script type="text/javascript">
    function toggle(obj) {
    var elstyle = document.getElementById(obj).style;
    var text    = document.getElementById(obj + "tog");
    if (elstyle.display == \'none\') {
    elstyle.display = \'block\';
    text.innerHTML = "hide";
    } else {
    elstyle.display = \'none\';
    text.innerHTML = "show";
    }
    }</script>
    </head>"""
    body = """<body>
    <p><a href="#Lods_Table" onclick="toggle(\'lods_tbl\')">Ranked lod scores</a>
    <div id="lods_tbl" class="divinfo">{}</div></p>
    <p><a href="#Hlods_Table" onclick="toggle(\'hlods_tbl\')">Ranked Hlod scores</a>
    <div id="hlods_tbl" class="divinfo">{}</div></p>
    <p><a href="#Lods_Heatmap" onclick="toggle(\'lods_heatmap\')">Lod scores heatmap</a>
    <div id="lods_heatmap">{}</div></p>
    <p><a href="#Hlods_Heatmap" onclick="toggle(\'hlods_heatmap\')">Hlod scores heatmap</a>
    <div id="hlods_heatmap">{}</div></p>
    </body>
    </html>"""
    with open('{0}/{0}_Report.html'.format(env.output), 'w') as f:
        #t = Template(index)
        #c = Context({ "lods": lods_tbl })
        f.write(head + body.format(html_table('Lod', theta_inc, theta_max, limit), html_table('Hlod', theta_inc, theta_max, limit), html_img('lod'), html_img('hlod')))

def html_img(ltype):
    chrs = ['{}'.format(i+1) for i in range(22)] + ['X', 'Y']
    imgs = ['heatmap/{0}.chr{1}.{2}s.png'.format(env.output, chrID, ltype) for chrID in chrs]
    exist_imgs = [img for img in imgs if isfile('{}/{}'.format(env.output,img))]
    return ''.join('<img src={}></img>'.format(img) for img in exist_imgs)
        
def html_table(type, theta_inc, theta_max, limit):
    colNum = int(round(theta_max/theta_inc))
    theta_orders = range(colNum)
    thetas = map(lambda x: theta_inc * x, range(colNum))
    #table
    table = r'<table style="width:300px;font-size:12px">{}</table>'
    #table header
    lods_header = r'<tr>{}</tr>'.format(''.join(r'<th colspan="2">&#952;:{}</td>'.format(x) for x in thetas))
    lods_header += r'<tr>{}</tr>'.format(r'<th rowspan="2">{}</th><th>Gene</th>'.format(type) * colNum)
    lods_header += r'<tr>{}</tr>'.format(r'<th>chr:start-end</th>' * colNum)
    #initialize lods dict
    lods = {}
    for theta_order in theta_orders:
        lods[theta_order] = {}
        #print '{}\n'.format(theta_order)
    #read lods
    lods_files = glob.glob('{0}/heatmap/{0}.*.{1}s'.format(env.output, type.lower()))
    for file in lods_files:
        with open(file, 'r') as f:
            for line in f:
                if type == 'Lod':
                    gene, chrId, start, end, theta, lod = line.strip().split()
                    if int(round(float(theta)/theta_inc)) >= colNum:
                        continue
                    lods[int(round(float(theta)/theta_inc))][gene] = [round(float(lod),3), gene, chrId, start, end]
                elif type == 'Hlod':
                    gene, chrId, start, end, a, theta, lod = line.strip().split()
                    if int(round(float(theta)/theta_inc)) >= colNum:
                        continue
                    lods[int(round(float(theta)/theta_inc))][gene] = ['{}<br>&#945;={}'.format(round(float(lod),3), round(float(a),3)), gene, chrId, start, end]
                else:
                    env.error('Wrong type of LOD')
                
    #collect result
    res = np.empty((len(theta_orders), limit), dtype=list)
    for theta_order in theta_orders:
        i=0
        for gene in sorted(lods[theta_order].keys(), key=lambda x: lods[theta_order][x][0], reverse=True):
            if i >= limit:
                break     
            res[theta_order][i] = lods[theta_order][gene]
            i += 1
    #print res
    #write lods table
    lods_res = ''
    for i in range(min(limit, len(lods[0].keys()))):
        lods_res += r'<tr>{}</tr>'.format(''.join(r'<td rowspan="2">{}</td><td>{}</td>'.format(res[theta_order][i][0],res[theta_order][i][1]) for theta_order in theta_orders))
        lods_res += r'<tr>{}</tr>'.format(''.join(r'<td>{}:{}-{}</td>'.format(res[theta_order][i][2],res[theta_order][i][3],res[theta_order][i][4]) for theta_order in theta_orders))
    lods_tbl = table.format(lods_header + lods_res)
    return lods_tbl
    
