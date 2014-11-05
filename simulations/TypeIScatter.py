#!/usr/bin/env python
import matplotlib.pyplot as plt

def plot(dfile, title):
    data = {'GJB2':[], 'SLC26A4':[]}
    with open(dfile) as f:
        for line in f.readlines():
            line = line.strip().split()
            data[line[0]].append(float(line[1]))
    plt.scatter([i+1 for i in range(len(data['SLC26A4']))], data['SLC26A4'])
    plt.xlim(0,500)
    plt.ylim(0,4.0)
    plt.axhline(y=3.6, color='r')
    plt.xlabel('Replicates')
    plt.ylabel('Cumulative HLOD score of 20 families')
    plt.title(title + '\n')

plt.figure(figsize=(12,12))
plt.subplot(3,1,1)
plot('TypeI/FAM20RUN0.hlods', '(A) Independent variants, no recombination, no missing data')
plt.subplot(3,1,2)
plot('TypeI/FAM20RUN1.hlods', '(B) Independent variants, within gene recombination, no missing data')
plt.subplot(3,1,3)
plot('TypeI/FAM20RUN0LD1PMissing.hlods', '(B) Variants complete LD, no recombination, 1 parent missing')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=2.0)
plt.savefig('typeI.png', dpi = 500)
