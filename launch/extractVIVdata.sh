cd /cluster/scratch_xp/public/cconti/CubismUP/Thesis/RisingVIV_varNu_15022016/
for f in $(ls *diagnostics.dat)
do
#	echo $f
	cat $f | awk '{ step=$1;t=$2;dt=$3;u=$10;v=$11;px+=dt*u;py+=dt*v;print step,t,px,py,u,v }' > $f.toPlot
done
cd -