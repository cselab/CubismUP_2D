#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 0.00001 -maxDT 0.00001
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 32 -maxBPD 32 -test advection -ic 1 -minDT 0.000001 -maxDT 0.0001 -tEnd 0.01

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=true -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 0.00001 -maxDT 0.00001
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 32 -maxBPD 32 -test advection -ic 1 -minDT 0.000001 -maxDT 0.0001 -tEnd 0.01

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=false particles=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 0.00001 -maxDT 0.00001
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 32 -maxBPD 32 -test advection -ic 1 -minDT 0.000001 -maxDT 0.0001 -tEnd 0.01

cd ../makefiles
make clean
make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=true particles=false -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 0.00001 -maxDT 0.00001
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 32 -maxBPD 32 -test advection -ic 1 -minDT 0.000001 -maxDT 0.0001 -tEnd 0.01