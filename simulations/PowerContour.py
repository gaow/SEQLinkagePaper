import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from collections import Counter
from tigerlib.utils import runCommand
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

def guess_seq_len(seq):
    guess = 1
    max_len = len(seq) / 2
    for x in range(2, max_len):
        if seq[0:x] == seq[x:2*x] :
            return x
    return guess

def to_matrix(l, n, transpose = False):
    out = [l[i:i+n] for i in range(0, len(l), n)]
    if transpose:
        return zip(*out)
    else:
        return out
    
class Plotter:
    def __init__(self, db):
        self.db = db
        self.moi = None
        self.allelic_het = 0
        self.id = 1
        self.name = None
        self.bands = 10
        self.colormap = {'SNV':"#0072B2", 'CHP':'#D55E00'}
        self.transparency = {'CHP':1, 'SNV':0.5}
        self.namemap = {'dominant':{1:'MYO7A', 2:'MYH9'}, 'recessive':{1:'SLC26A4', 2:'GJB2'},
                        'compound_recessive':{1:'SLC26A4', 2:'GJB2'}}
        self.thresholdmap = {'lod':3.3, 'hlod':3.6}
        
    def SetParam(self, i, moi, a, score):
        self.id = i
        self.moi = moi
        self.allelic_het = a
        self.name = self.namemap[self.moi][self.id]
        self.score = score

    def GetData(self, tables, group_by = None):
        data = {}
        for table in tables:
            where = 'where moi = \'{}\''.format(self.moi) if self.moi else ''
            where += (" AND " if where else '') + 'ahet = \'{}\''.format(self.allelic_het)
            cmd = 'wsqlite {0} "select fam_size, prop{1}, p{4}{1} from {2} {3} order by fam_size, prop{1}"'.\
              format(self.db, self.id, table, where, self.score)
            out = zip(*[map(float, x.split()) for x in runCommand(cmd).strip().split('\n')])
            if group_by is None:
                data[table] = out
            else:
                values = Counter(out[group_by - 1]).values()
                if len(set(values)) != 1:
                    raise ValueError('Matrix rows not equal')
                # unit = len(out[group_by - 1]) / values[0]
                data[table] = [to_matrix(x, values[0], transpose = True) for x in out]
        return data

    def Plot(self, data, out = 'contour.pdf'):
        plt.figure(figsize=(10,10))
        for t in data:
            cs = plt.contour(data[t][0], data[t][1], data[t][2], self.bands + 1, colors = self.colormap[t],
                             alpha = self.transparency[t])
            plt.clabel(cs, inline=1, fontsize=12, fmt='%1.2f')
        plt.title("{} [{}{}{}]\n".format(self.name, self.score.upper(), r'$\geq$',
                                         self.thresholdmap[self.score]), fontsize = 20)
        plt.xlabel("Number of families", fontsize = 20)
        plt.ylabel("Locus heterogeneity\n", fontsize = 20)
        plt.savefig(out, dpi = 500)


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        sys.argv.append('lod')
        task = 'lod'
    else:
        task = sys.argv[1]
    p = Plotter('PowerCalc.sqlite3')
    for moi in ['recessive', 'dominant', 'compound_recessive']:
        for i in [1,2]:
            for a in [0, 1]:
                p.SetParam(i, moi, a, sys.argv[-1])
                p.Plot(p.GetData(['CHP', 'SNV'], group_by = 1), out = "PowerFigs/{}_a{}_{}_gene{}.pdf".format(task,a,moi,i))
