import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

polygon = [(0.0, 0.0), (1.0, 1.0), (2.0, 1.0),
           (3.0, 0.5), (3.0, -1.0), (2.0, -2.5),
           (1.0, -2.5), (0.0, -1.5), (1.0, -0.5)]

def Ruppert(points, segments, threshold):
    pass

def TriangPolygon(p):
    # Simple, but not optimal triangulation of a polygon
    segments = []
    points = p
    n = len(polygon)
    for i in range(0, n-1):
        segments.append([p[i], p[i+1]])
    segments.append([p[n-1] ,p[0]])
    for i in range(0, int((n-1)/2)): # "Downwards" segments
        segments.append([p[i+1], p[n-1-i]])
    #for i in range(0, int(n/2)-1):
    #    segments.append([p[n-1-i], p[i+2]]) # "Upwards" segments
    return points, segments

def FlipDelauney(points, segments):
    # Delauney triangulation based on the flip algorithm
    # O(n^2)
    return points, segments

def DrawTriangulation(points, segments):
    plt.figure(0)
    ax = plt.gca()
    lc = LineCollection(segments, linewidths=1, linestyle="solid")
    ax.add_collection(lc)
    ax.set_xlim(-0.5, 4.0)
    ax.set_ylim(-3.5, 1.5)
    plt.show()
    
def main():
    p, s = TriangPolygon(polygon)
    print(s)
    DrawTriangulation(p, s)
    
if __name__ == '__main__':
    main()