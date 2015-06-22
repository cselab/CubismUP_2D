cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed precision=double particles=false -j
cd ../launch
export OMP_NUM_THREADS=48

for R in 0.01 0.1 0.5 0.9 1.0 1.1 2.0 10.0 100.0 1000.0 10000.0
do
	../makefiles/test -file /cluster/scratch_xp/public/cconti/CubismUP/Disk${I}_rhoS${R} -CFL 0.1 -LCFL 0.1 -minBPD 4 -maxBPD 64 -test addedmass -rhoS ${R} -radius 0.05 -nu 0.001
#	../makefiles/test -file /cluster/scratch_xp/public/cconti/CubismUP/Disk_split${I}_rhoS${R} -CFL 0.1 -LCFL 0.1 -minBPD 4 -maxBPD 64 -test addedmass -rhoS ${R} -radius 0.05 -nu 0.001 -split
done