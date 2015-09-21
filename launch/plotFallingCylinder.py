import os
import math
import numpy as np
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib

fig = figure()
fig.canvas.set_window_title('Falling Cylinders')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'

data = []
x = []
y = []
r = []
v = []
xref = [0, 20]
rref = [40, 40]
vref = [-.48, -.48]

p1 = fig.add_subplot(131)
p2 = fig.add_subplot(132)
p3 = fig.add_subplot(133)

for dirName, subDirList, fileList in os.walk(rootDir):
	if "Falling_Re40_0409" in dirName:
		for file in fileList:
			if "diagnostics.dat" in file:
				#and "128" in file:
				data.append(np.genfromtxt(fname=dirName+'/'+file))
				idx = len(data)-1
				dataset = data[idx]
				x.append(dataset[:,1])
				r.append(dataset[:,6])
				y.append(dataset[:,8])
				v.append(dataset[:,10])
				p1.plot(x[idx],r[idx],label=file)
				p2.plot(x[idx],v[idx],label=file)
				p3.plot(x[idx],y[idx],label=file)

				p3.legend()

p1.plot(xref,rref)
p2.plot(xref,vref)

handles, labels = p3.get_legend_handles_labels()
p3.legend(handles, labels)

p1.set_xlim(0,2)
p1.set_ylim(0,100)
p2.set_xlim(0,2)
p2.set_ylim(-1,0)
p3.set_xlim(0,2)
p3.set_ylim(0,1)

p1.set_title("Re40")
p1.set_xlabel('Time')
p1.set_ylabel('Re')

p2.set_title("Re40")
p2.set_xlabel('Time')
p2.set_ylabel('velocity v')

p3.set_title("Re40")
p3.set_xlabel('Time')
p3.set_ylabel('position y')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()