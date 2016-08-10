awk '{ step=$1;t=$2;theta=$12;omega=$13;if (theta>atan2(0, -1)) theta=theta-2*atan2(0, -1);print step,t,theta,omega }' VarDensityDisk_TopLight_1102/VarDensityDisk_TopLight_1102_diagnostics.dat
