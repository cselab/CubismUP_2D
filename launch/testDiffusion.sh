cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test diffusion