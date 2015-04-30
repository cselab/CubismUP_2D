cd ../makefiles/
make clean
make bc=periodic config=production multiphase=false precision=single vertexcentered=true poisson=fftw -j
cd -
export OMP_NUM_THREADS=48
../makefiles/test -minBPD 4 -maxBPD 4 -test shearlayer
#../makefiles/test -minBPD 1 -maxBPD 16 -test shearlayer
