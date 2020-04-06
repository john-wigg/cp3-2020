import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

file = open("build-lin/delaunay.out", "r")
lines = file.readlines()
print(lines)#
line = lines[0]
i = 1

vertices = []
triangles = []
segments = []

while line != '\n':
    vertices.append(list(map(float, line.rstrip().split(" "))))
    line = lines[i]
    i+=1

line = lines[i]  
i+=1  
while line != '\n':
    triangles.append(list(map(int, line.rstrip().split(" "))))
    line = lines[i]
    i+=1
    
line = lines[i]  
i+=1  
while line != '\n':
    segments.append(list(map(int, line.rstrip().split(" "))))
    line = lines[i]
    i+=1
    
print(vertices)
print(triangles)
print(segments)

plt.figure(0, figsize=(10, 10))
ax = plt.gca()

segs = []

for i in range(len(triangles)):
    x = 0
    y = 0
    for k in range(3):
            line = (vertices[triangles[i][k]], vertices[triangles[i][(k+1) % 3]])
            x += vertices[triangles[i][k]][0]
            y += vertices[triangles[i][k]][1]
            segs.append(line)
    x /= 3
    y /= 3
    plt.text(x, y, str(i), color='grey')
    
lc = LineCollection(segs, linewidths=1, linestyle="solid")
ax.add_collection(lc)

for i in range(len(vertices)):
    x = vertices[i][0]
    y = vertices[i][1]
    plt.text(x, y, str(i))

segs_poly = []
for s in segments:
    line = (vertices[s[0]], vertices[s[1]])
    segs_poly.append(line)
    
lc_poly = LineCollection(segs_poly, linewidths=1, color='red')
ax.add_collection(lc_poly)

ax.set_xlim(-5, 5)
ax.set_ylim(-6.5, 5)
ax.set_aspect('equal')

for vertex in vertices:
    plt.plot(vertex[0], vertex[1], "x")

#c = plt.Circle((3.42402, -4.69935), 5.39287)
#ax.add_artist(c)

plt.show()