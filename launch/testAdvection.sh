#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 multiphase=true -j
#make config=debug poisson=hypre bc=vortex precision=double bs=32 multiphase=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 32 -test advection -ic 0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 2 -test advection -ic 0

cd ../makefiles
make clean
make config=production poisson=hypre bc=vortex precision=double bs=32 multiphase=false -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 4 -test advection -ic 1
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 32 -test advection -ic 1