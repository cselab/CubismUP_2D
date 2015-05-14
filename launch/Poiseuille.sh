cd ../makefiles/
make clean
make bc=pipe config=production multiphase=false precision=single poisson=fftw particles=true -j
cd -
export OMP_NUM_THREADS=48
mpirun -np 1 ../makefiles/test -minBPD 4 -maxBPD 4 -test poiseuille -file ../data/testPoiseuille

