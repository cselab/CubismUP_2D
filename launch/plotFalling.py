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
xref = [0, 2]
vrefk = [-.53, -.53]
vrefm = [-.47, -.47]

pku = fig.add_subplot(221)
pkp = fig.add_subplot(222)
pmu = fig.add_subplot(223)
pmp = fig.add_subplot(224)

for dirName, subDirList, fileList in os.walk(rootDir):
#if "Kolomenskiy" in dirName:
	for file in fileList:
		if "diagnostics.dat" in file and "Kolomenskiy" in file:
			data.append(np.genfromtxt(fname=dirName+'/'+file))
			idx = len(data)-1
			dataset = data[idx]
			x.append(dataset[:,1])
			r.append(dataset[:,6])
			y.append(dataset[:,8])
			v.append(dataset[:,10])
			pku.plot(x[idx],v[idx],label=file)
			pkp.plot(x[idx],y[idx],label=file)
#if "Mordant" in dirName:
	for file in fileList:
		if "diagnostics.dat" in file and "Mordant" in file:
			data.append(np.genfromtxt(fname=dirName+'/'+file))
			idx = len(data)-1
			dataset = data[idx]
			x.append(dataset[:,1])
			r.append(dataset[:,6])
			y.append(dataset[:,2])
			#y.append(dataset[:,8])
			v.append(dataset[:,10])
			pmu.plot(x[idx],v[idx],label=file)
			pmp.plot(x[idx],y[idx],label=file)
			pmu.legend()


handles, labels = pmu.get_legend_handles_labels()
pmu.legend(handles, labels)

pku.plot(xref,vrefk)
pmu.plot(xref,vrefm)

pku.set_xlim(0,2)
pku.set_ylim(-1,0)
pkp.set_xlim(0,2)
pkp.set_ylim(0,1)

pmu.set_xlim(0,2)
pmu.set_ylim(-1,0)
pmp.set_xlim(0,2)
pmp.set_ylim(0,0.001)
#pmp.set_ylim(0,1)

pku.set_title("Kolomenskiy")
pku.set_xlabel('Time')
pku.set_ylabel('velocity v')

pkp.set_title("Kolomenskiy")
pkp.set_xlabel('Time')
pkp.set_ylabel('position y')

pmu.set_title("Mordant")
pmu.set_xlabel('Time')
pmu.set_ylabel('velocity v')

pmp.set_title("Mordant")
pmp.set_xlabel('Time')
pmp.set_ylabel('position y')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()