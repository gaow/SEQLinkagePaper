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
    out = [l[i:i+n] for i in xrange(0, len(l), n)]
    if transpose:
        return zip(*out)
    else:
        return out
    
class Plotter:
    def __init__(self, db):
        self.db = db
        self.moi = None
        self.id = 1
        self.bands = 10
        self.colormap = {'SNV':"#0072B2", 'CHP':'#D55E00'}
        self.transparency = {'CHP':1, 'SNV':0.5}
        
    def SetMOI(self, moi):
        self.moi = moi

    def SetID(self, i):
        self.id = i

    def GetData(self, tables, group_by = None):
        data = {}
        for table in tables:
            where = 'where moi = \'{}\''.format(self.moi) if self.moi else ''
            cmd = 'wsqlite {0} "select fam_size, prop{1}, plod{1} from {2} {3} order by fam_size, prop{1}"'.\
              format(self.db, self.id, table, where)
            out = zip(*[map(float, x.split()) for x in runCommand(cmd).strip().split('\n')])
            if group_by is None:
                data[table] = out
            else:
                values = Counter(out[group_by - 1]).values()
                if len(set(values)) != 1:
                    raise ValueError('Matrix rows not equal')
                unit = len(out[group_by - 1]) / values[0]
                data[table] = [to_matrix(x, unit, transpose = True) for x in out]
        return data

    def Plot(self, data, out = 'contour.pdf'):
        plt.figure(figsize=(10,10))
        for t in data:
            cs = plt.contour(data[t][0], data[t][1], data[t][2], self.bands + 1, colors = self.colormap[t],
                             alpha = self.transparency[t])
            plt.clabel(cs, inline=1, fontsize=12, fmt='%1.2f')
        plt.title("Gene{}, {}\n".format(self.id, self.moi), fontsize = 20)
        plt.xlabel("Family size", fontsize = 20)
        plt.ylabel("Heterogeneity\n", fontsize = 20)
        plt.savefig(out, dpi = 500)


if __name__ == '__main__':
    p = Plotter('PowerCalc.sqlite3')
    for moi in ['recessive', 'dominant', 'compound_recessive', 'compound_dominant']:
        for i in [1,2]:
            p.SetMOI(moi)
            p.SetID(i)
            p.Plot(p.GetData(['CHP', 'SNV'], group_by = 1), out = "PowerFigs/{}_gene{}.pdf".format(moi,i))
