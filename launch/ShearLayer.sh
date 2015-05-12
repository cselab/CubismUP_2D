cd ../makefiles/
make clean
make bc=periodic config=production multiphase=false precision=single vertexcentered=true poisson=fftw -j
cd -
export OMP_NUM_THREADS=48
../makefiles/test -minBPD 8 -maxBPD 8 -test shearlayer -file ../data/testShearLayer_Euler_Mp4_
#../makefiles/test -minBPD 1 -maxBPD 16 -test shearlayer -file ../data/testShearLayer_Convergence_
