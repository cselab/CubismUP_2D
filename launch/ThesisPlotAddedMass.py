import os
import math
import numpy as np
from pylab import figure, show
import matplotlib.pyplot as plt
import matplotlib

fig = figure()
fig.canvas.set_window_title('Falling Cylinders')

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/Thesis/'

data = []
x = []
y = []

p = fig.add_subplot(111)

for dirName, subDirList, fileList in os.walk(rootDir):
	if "AddedMass" in dirName:
		for file in fileList:
			if "addedmass.dat" in file:
				data.append(np.genfromtxt(fname=dirName+'/'+file))
				idx = len(data)-1
				dataset = data[idx]
				x.append(dataset[:,1])
				y.append(dataset[:,5])
				p.plot(x[idx],y[idx],label=file)


#handles, labels = p3.get_legend_handles_labels()
#p3.legend(handles, labels)


p.set_title("Added Mass")
p.set_xlabel('dh')
p.set_ylabel('Relative Error')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
show()