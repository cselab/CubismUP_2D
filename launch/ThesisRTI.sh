for p in true
#false#true
do
for rk in true false
do
cd ../makefiles/
make clean
make bc=mixed config=production multiphase=true precision=double vertexcentered=true poisson=split-fftw particles=${p} densitydiff=true vardensity=true rk2=${rk} -j
cd -

BASENAME=Thesis_RTI_25032016_Particles${p}_RK2${rk}_variableRho
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
mkdir ${FOLDER}

for BPD in 16 32
#2 4 8
do
cp ../makefiles/simulation ${FOLDER}/simulation_BPD${BPD}_Particles${P}_RK2${rk}
cd ${FOLDER}
export OMP_NUM_THREADS=48
bsub -n 48 -W 10:00 -o ${BASENAME}_BPD${BPD}_Particles${P}_RK2${rk} -J ${BASENAME}_BPD${BPD}_Particles${P}_RK2${rk} ./simulation_BPD${BPD}_Particles${P}_RK2${rk} -file RTI_BPD${BPD} -CFL 0.1 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx ${BPD} -bpdy $[4*${BPD}] -tend 1. -sim rti -fdump $[20*${BPD}] -nsteps $[2000*${BPD}]
cd -
done
done
done
