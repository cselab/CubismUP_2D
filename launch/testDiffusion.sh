cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/

export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -minDT 1e-8 -maxDT 1e-8 -test diffusion
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -minDT 1e-9 -maxDT 1e-7 -test diffusion