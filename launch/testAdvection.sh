#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 multiphase=true -j
#make config=debug poisson=hypre bc=vortex precision=double bs=32 multiphase=true -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 32 -test advection -ic 0 -minDT 0.000001 -maxDT 0.000001
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 0 -minDT 0.000001 -maxDT 0.0001 -tEnd 0.01

cd ../makefiles
make clean
make config=production poisson=hypre bc=vortex precision=double bs=32 -j
#make config=debug poisson=hypre bc=vortex precision=double bs=32 -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 0.00001 -maxDT 0.00001
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 1 -minDT 0.000001 -maxDT 0.0001 -tEnd 0.01