cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/

export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -minDT 1e-8 -maxDT 1e-8 -test diffusion
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-6 -maxDT 1e-4 -test diffusion