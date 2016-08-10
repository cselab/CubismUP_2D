cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed precision=double particles=false movingframe=false dlm=true -j
cd ../launch
export OMP_NUM_THREADS=48

BASENAME=AddedMass_DP_18022016_inOut_rhoS
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
mkdir ${FOLDER}

cp ../makefiles/test ${FOLDER}/
cd ${FOLDER}
for R in 0.01 0.1 0.5 0.9 1.0 1.1 2.0 10.0 100.0 1000.0 10000.0
do
	bsub -n 48 -W 1:00 -o ${BASENAME}_Rho${R} -J ${BASENAME}_Rho${R} ./test -file AddedMass_Cylinder_rhoS${R} -CFL 0.1 -LCFL 0.1 -minBPD 4 -maxBPD 64 -test addedmass -rhoS ${R} -radius 0.05 -nu 0.001 -nsteps 1 dlm 1
done