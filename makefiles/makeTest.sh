make clean;make -j
make clean;make bc=periodic -j
make clean;make precision=double -j
make clean;make poisson=split-fftw
make clean;make multiphase=true

make clean;make config=production -j
make clean;make config=production bc=periodic -j
make clean;make config=production precision=double -j
make clean;make config=production poisson=split-fftw
make clean;make config=production multiphase=true