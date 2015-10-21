cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=double particles=false rk2=true -j
cd ../launch

for BPD in 16 32 64
do
	BASENAME=Andersen_2110_CFL0.01_diff4th_bpd
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}${BPD}
    mkdir ${FOLDER}
    cp ../makefiles/simulation ${FOLDER}/simulation_T
    cp ../makefiles/simulation ${FOLDER}/simulation_F
    cp Andersen.sh ${FOLDER}
    cd ${FOLDER}
    export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}${BPD} -J ${BASENAME}${BPD} ./simulation_T -file ${FOLDER}/${BASENAME}${BPD}_T -CFL 0.01 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 2.7 -ypos .9 -xpos .5 -nu 0.0000190757887517147 -sim falling -shape ellipse -semiAxisX 0.003125 -semiAxisY 0.025 -angle 0.200712863979348 -tdump 0.1 -lambda 1e6
    export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}${BPD} -J ${BASENAME}${BPD} ./simulation_F -file ${FOLDER}/${BASENAME}${BPD}_F -CFL 0.01 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.35 -ypos .9 -nu 0.0000190757887517147 -sim falling -shape ellipse -semiAxisX 0.003125 -semiAxisY 0.025 -angle 0.200712863979348 -tdump 0.1 -lambda 1e6
    cd -
done