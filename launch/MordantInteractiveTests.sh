cd ../makefiles
make clean
make config=debug poisson=hypre bc=mixed precision=double particles=false -j
cd ../launch

BPD=16

FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Mordant_test2_bpd${BPD}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
export OMP_NUM_THREADS=48;mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/Mordant_test2_bpd${BPD}/Mordant_test2_bpd${BPD} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -radius .01 -tend 2. -rhoS 2.565 -ypos .9 -nu 0.000228609536213148 -sim falling -tdump 0.01 -lambda 1e4
