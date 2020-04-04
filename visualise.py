import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

pos = np.genfromtxt("square.in", skip_header=2, skip_footer=9)
disp = np.loadtxt("square.out")

num_nodes = int(disp.size / 2)

fig = plt.figure(0)
ax = plt.gca()

segs = []
segs_no_disp = []

for i in range(num_nodes):
    for k in range(num_nodes):
        line = (disp[2*i] + pos[i,0], disp[2*i+1] + pos[i, 1]), (disp[2*k] + pos[k, 0], disp[2*k+1] + pos[k, 1])
        segs.append(line)
        line = (pos[i,0], pos[i, 1]), (pos[k, 0], pos[k, 1])
        segs_no_disp.append(line)
    
lc = LineCollection(segs, linewidths=1, linestyle="solid")
ax.add_collection(lc)

lc = LineCollection(segs_no_disp, linewidths=0.5, color="black", linestyle="dashed")
ax.add_collection(lc)

ax.set_xlim(-0.2, 1.2)
ax.set_ylim(-0.2, 1.2)


plt.show()