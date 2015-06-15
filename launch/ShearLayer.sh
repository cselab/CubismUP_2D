cd ../makefiles/
make clean
make bc=periodic config=production multiphase=false precision=single vertexcentered=true poisson=fftw particles=false -j
cd -
export OMP_NUM_THREADS=48
../makefiles/test -minBPD 4 -maxBPD 4 -test shearlayer -file ../data/testShearLayer_Euler_FD_
#../makefiles/test -minBPD 1 -maxBPD 16 -test shearlayer -file ../data/testShearLayer_Convergence_
