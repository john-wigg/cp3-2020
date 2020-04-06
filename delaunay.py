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
    
print(vertices)
print(triangles)

plt.figure(0)
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
    plt.text(x, y, str(i))
    
lc = LineCollection(segs, linewidths=1, linestyle="solid")
ax.add_collection(lc)

ax.set_xlim(-3.5, 3.5)
ax.set_ylim(-3.5, 3.5)
ax.set_aspect('equal')

for vertex in vertices:
    plt.plot(vertex[0], vertex[1], "x")

#c = plt.Circle((3.42402, -4.69935), 5.39287)
#ax.add_artist(c)

plt.show()