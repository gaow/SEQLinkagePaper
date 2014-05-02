import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import prettyplotlib as ppl
import numpy as np
import string

fig, ax = plt.subplots(1)

np.random.seed(10)

ppl.pcolormesh(fig, ax, np.random.randn(10,10), 
				xticklabels=np.array(range(10)), 
		        yticklabels=string.lowercase[-10:])
fig.savefig('pcolormesh_prettyplotlib_labels.png')
