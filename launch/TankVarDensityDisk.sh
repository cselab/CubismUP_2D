cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=double rk2=true particles=false dlm=true constnu=true tanktreading=true -j
#
cd ../launch

NAME=Tanktreading_VarDensityDisk_.2PIOscillatingOmega_1512
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp TankVarDensityDisk.sh ${FOLDER}
cd ${FOLDER}
export OMP_NUM_THREADS=48
#bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME}
./simulation -file ${FOLDER}/${NAME} -CFL 0.01 -dlm 1 -lambda 1e6 -bpdx 16 -bpdy 16 -radius .1 -tend 1000. -rhoS 0.99 -rhoS2 1.01 -angle .0 -ypos .5 -xpos .5 -nu 0.0001 -sim falling -shape diskVarDensity -fdump 100 -eps 2
cd -