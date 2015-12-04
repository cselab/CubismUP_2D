import os
import math
import numpy as np
from matplotlib.patches import Ellipse
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec

fig = figure(figsize=(16,9))
fig.canvas.set_window_title('Falling Ellipses')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'

data = []
u = []
v = []
t = []
dtdt = []
x = []
y = []
a = []
n = []
ells = []

gs = gridspec.GridSpec(2,3)
p  = fig.add_subplot(gs[:,:-1],aspect='equal')
pu = fig.add_subplot(gs[0,2],aspect='equal')
pa = fig.add_subplot(gs[1,2])#,aspect='equal')
colorList = ['#0088d9','#d95100','#00d951','#d90088']
counter = 0

for dirName, subDirList, fileList in os.walk(rootDir):
	for file in fileList:
		if "Andersen_0412_Re1100_small_box_DLM0.75_CFL0.01_bpd16_.1_T_diagnostics.dat" in file and "keep" in dirName:
			data.append(np.genfromtxt(fname=dirName+'/'+file))
			idx = len(data)-1
			dataset = data[idx]
			u.append(dataset[:,9])
			v.append(dataset[:,10])
			t.append(dataset[:,11])
			dtdt.append(dataset[:,12])
			pu.plot(u[idx],v[idx],label=file,color=colorList[counter],linewidth=2.0)
			pa.plot(t[idx]%(2.*math.pi),dtdt[idx],label=file,color=colorList[counter],linewidth=0.0,marker='.')
				
			x.append(dataset[:,7])
			y.append(dataset[:,8])
			a.append(dataset[:,11])
			n.append(x[idx].size)
			#ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.2, height=0.025, angle=a[idx][i]*360/(2*math.pi))
			#ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.1, height=0.0125, angle=a[idx][i]*360/(2*math.pi))
			#ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.05, height=0.00625, angle=a[idx][i]*360/(2*math.pi))
			ells.append([Ellipse(xy=(x[idx][i],y[idx][i]), width=0.025, height=0.003125, angle=a[idx][i]*360/(2*math.pi))
						 for i in range(n[idx])])
			increment = 10
			for i in range(1,n[idx],increment):
				e = ells[idx][i]
				p.add_artist(e)
				#p.legend()
				e.set_clip_box(p.bbox)
				e.set_alpha(float(i)/float(n[idx]))
				e.set_facecolor(colorList[counter])
#e.set_facecolor((float(bIdx)*.25,float(3-bIdx)*.05,float(3-bIdx)*.25))
			counter=(counter+1)%4

pu.set_xlim(-.75,.75)
pu.set_ylim(-.75,.75)
pa.set_xlim(0,2*math.pi)
pa.set_ylim(-20,20)
#pa.set_xlim(-10,10)
#pa.set_ylim(-10,10)

pu.set_xlabel('u')
pu.set_ylabel('v')
pa.set_xlabel('theta')
pa.set_ylabel('dthetadt')

p.set_xlim(0, 1)
p.set_ylim(0, 1)
p.set_xlabel('x')
p.set_ylabel('y')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
fig.savefig("Andersen.png",dpi=300,transparent=True)
show()