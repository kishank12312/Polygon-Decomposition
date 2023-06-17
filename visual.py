#DAA assignment visualizer
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)
def display(filename):
    f = open(filename)
    numpoly = int(f.readline())
    polygons = []
    for i in range(numpoly):
        f.readline()
        polysize = int(f.readline())
        polygon = []
        for j in range(polysize):
            x,y = map(float,f.readline().split())
            polygon.append((x,y))
        polygons.append(polygon)
    minx,miny,maxx,maxy = 1000000000000000,1000000000000000,-1000000000000000,-1000000000000000
    for i in polygons:
        for j in i:
            x,y = j
            minx = min(minx,x)
            miny = min(miny,y)
            maxx = max(maxx,x)
            maxy = max(maxy,y)
    cmap = get_cmap(len(polygons))
    fig, ax = plt.subplots()
    for i in polygons:
        for j in i:
            x,y = j
            plt.plot(x,y, marker = 'o', color = 'r')

    for i,poly in enumerate(polygons):
        ax.add_patch(Polygon(poly, facecolor=cmap(i), edgecolor = 'k',linestyle='-'))
        x = sum([i[0] for i in poly])/len(poly)
        y = sum([i[1] for i in poly])/len(poly)
        ax.text(x,y, str(i), ha ='center',va ='center', size = 20)

    ax.set_xlim([minx-(0.1)*(maxx-minx), maxx+(0.1)*(maxx-minx)])
    ax.set_ylim([miny-(0.1)*(maxy-miny), maxy+(0.1)*(maxy-miny)])

    plt.show()

display("Unmerged Polygon.txt")
display("Merged Polygon.txt")