for p in false
#true
do
for rk in false
#true
do
cd ../makefiles/
make clean
make bc=periodic config=production multiphase=false precision=double poisson=fftw rk2=${rk} particles=${p} -j
cd -

BASENAME=Thesis_TravelingWave_220316_Particles${p}_RK2${rk}
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
mkdir ${FOLDER}

cp ../makefiles/test ${FOLDER}/test_Particles${P}_RK2${rk}
cd ${FOLDER}
export OMP_NUM_THREADS=48
bsub -n 48 -W 08:00 -o ${BASENAME} -J ${BASENAME} ./test_Particles${P}_RK2${rk} -minBPD 1 -maxBPD 16 -minDT 1 -maxDT 1 -test travelingwave -file TravelingWave
cd -
done
done
