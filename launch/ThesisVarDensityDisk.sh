cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=double rk2=true particles=false dlm=true constnu=true -j
cd ../launch

NAME=VarDensityDisk_TopLight_2002
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${NAME}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp ThesisVarDensityDisk.sh ${FOLDER}
cd ${FOLDER}
export OMP_NUM_THREADS=48
bsub -n 48 -W 40:00 -o ${NAME} -J ${NAME} ./simulation -file ${FOLDER}/${NAME} -CFL 0.01 -dlm 1 -bpdx 32 -bpdy 32 -radius .025 -tend 10. -rhoS 1.1 -rhoS2 1.2 -angle .79 -ypos .9 -xpos .5 -nu 0.00001 -sim falling -shape diskVarDensity -tdump 0.1 -eps 2
#bsub -n 48 -W 40:00 -o ${NAME} -J ${NAME} ./simulation -file ${FOLDER}/${NAME} -CFL 0.01 -dlm 1 -lambda 1e6 -bpdx 32 -bpdy 32 -radius .025 -tend 10. -rhoS 1.375 -rhoS2 .75 -angle .3 -ypos .5 -xpos .5 -nu 0.0000001 -sim falling -shape fish -tdump .05 -eps 2
cd -

NAME=VarDensityDisk_TopHeavy_2002
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${NAME}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp ThesisVarDensityDisk.sh ${FOLDER}
cd ${FOLDER}
export OMP_NUM_THREADS=48
bsub -n 48 -W 40:00 -o ${NAME} -J ${NAME} ./simulation -file ${FOLDER}/${NAME} -CFL 0.01 -dlm 1 -bpdx 32 -bpdy 32 -radius .025 -tend 10. -rhoS 1.2 -rhoS2 1.1 -angle .79 -ypos .9 -xpos .5 -nu 0.00001 -sim falling -shape diskVarDensity -tdump 0.1 -eps 2
cd -
