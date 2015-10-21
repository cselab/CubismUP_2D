import os
import math
import numpy as np
from matplotlib.patches import Ellipse
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 4})

fig = figure(dpi=280)
fig.canvas.set_window_title('Falling Ellipses (Andersen)')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'

data = []
x = []
y = []
a = []
n = []
ells = []

p = fig.add_subplot(111, aspect='equal')

bIdx = 0
for B in [32, 64]:
	dirName = 'Andersen_2110_CFL0.01_bpd'+str(B)
    #	dirName = 'Andersen_2110_bpd'+str(B)
	fileNameT = dirName+'_T_diagnostics.dat'
	fullNameT = rootDir+dirName+'/'+fileNameT
	if os.path.isfile(fullNameT):
		data.append(np.genfromtxt(fname=fullNameT))
		idx = len(data)-1
		dataset = data[idx]
		x.append(dataset[:,7])
		y.append(dataset[:,8])
		a.append(dataset[:,11])
		n.append(x[idx].size)
		ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.05, height=0.00625, angle=a[idx][i]*360/(2*math.pi))
					 for i in range(n[idx])])
		increment = 10
		for i in range(1,n[idx],increment):
			e = ells[idx][i]
			p.add_artist(e)
			e.set_clip_box(p.bbox)
			e.set_alpha(float(i)/float(n[idx]))
			e.set_facecolor((float(bIdx)*.33,float(3-bIdx)*.05,float(3-bIdx)*.33))

	fileNameF = dirName+'_F_diagnostics.dat'
	fullNameF = rootDir+dirName+'/'+fileNameF
	if os.path.isfile(fullNameF):
		data.append(np.genfromtxt(fname=fullNameF))
		idx = len(data)-1
		dataset = data[idx]
		x.append(dataset[:,7])
		y.append(dataset[:,8])
		a.append(dataset[:,11])
		n.append(x[idx].size)
		ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.05, height=0.00625, angle=a[idx][i]*360/(2*math.pi))
					 for i in range(n[idx])])
		increment = 10
		for i in range(1,n[idx],increment):
			e = ells[idx][i]
			p.add_artist(e)
			e.set_clip_box(p.bbox)
			e.set_alpha(float(i)/float(n[idx]))
			e.set_facecolor((float(bIdx)*.33,float(3-bIdx)*.25,float(3-bIdx)*.33))

	bIdx = bIdx+1

p.set_xlim(0, 1)
p.set_ylim(0, 1)

p.set_xlabel('x')
p.set_ylabel('y')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()

