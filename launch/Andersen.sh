cd ../makefiles
make clean
make config=production poisson=split-fftw bc=box precision=double particles=false rk2=true dlm=true -j
cd ../launch

for D in 0.5 0.75 1.0 1.25 1.5 1.75 2.0 5.0
do
for P in .1 .09 .11 .15 .2
do
for CFL in 0.01
#0.05 0.1 0.2
do
for BPD in 16 32
#32 64 128
do
	BASENAME=Andersen_0412_Re1100_box_DLM${D}_CFL${CFL}_bpd${BPD}_${P}
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
    mkdir ${FOLDER}
    cp ../makefiles/simulation ${FOLDER}/simulation_T
#    cp ../makefiles/simulation ${FOLDER}/simulation_F
    cp Andersen.sh ${FOLDER}
	cd ${FOLDER}
#export OMP_NUM_THREADS=48;bsub -n 48 -W 40:00 -o ${BASENAME}_T -J ${BASENAME}_T ./simulation_T -file ${FOLDER}/${BASENAME}_T -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 20. -rhoS 4.0 -ypos .9 -xpos ${P} -nu 0.0000190757887517147 -sim falling -shape ellipse -semiAxisX 0.003125 -semiAxisY 0.025 -angle 0.2 -tdump 0.02 -dlm 1
#export OMP_NUM_THREADS=48;bsub -n 48 -W 40:00 -o ${BASENAME}_F -J ${BASENAME}_F ./simulation_F -file ${FOLDER}/${BASENAME}_F -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 20. -rhoS 1.25 -ypos .9 -xpos .5 -nu 0.0000190757887517147 -sim falling -shape ellipse -semiAxisX 0.003125 -semiAxisY 0.025 -angle 0.2 -tdump 0.02 -dlm 1
#export OMP_NUM_THREADS=48;bsub -n 48 -W 40:00 -o ${BASENAME}${BPD}_T -J ${BASENAME}${BPD}_T ./simulation_T -file ${FOLDER}/${BASENAME}${BPD}_T -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 2.7 -ypos .9 -xpos .1 -nu 0.000152606310013717 -sim falling -shape ellipse -semiAxisX 0.0125 -semiAxisY 0.1 -angle 0.2 -tdump 0.1 -lambda 1e6
#export OMP_NUM_THREADS=48;bsub -n 48 -W 40:00 -o ${BASENAME}${BPD}_F -J ${BASENAME}${BPD}_F ./simulation_F -file ${FOLDER}/${BASENAME}${BPD}_F -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.35 -ypos .9 -nu 0.000152606310013717 -sim falling -shape ellipse -semiAxisX 0.0125 -semiAxisY 0.1 -angle 0.2 -tdump 0.1 -lambda 1e6
#export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}${BPD}_T -J ${BASENAME}${BPD}_T ./simulation_T -file ${FOLDER}/${BASENAME}${BPD}_T -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 2.7 -ypos .9 -xpos .1 -nu 0.0000539544783312781 -sim falling -shape ellipse -semiAxisX 0.00625 -semiAxisY 0.05 -angle 0.2 -tdump 0.1 -lambda 1e6
#export OMP_NUM_THREADS=48;bsub -n 48 -W 40:00 -o ${BASENAME}${BPD}_F -J ${BASENAME}${BPD}_F ./simulation_F -file ${FOLDER}/${BASENAME}${BPD}_F -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.35 -ypos .9 -nu 0.0000539544783312781 -sim falling -shape ellipse -semiAxisX 0.00625 -semiAxisY 0.05 -angle 0.2 -tdump 0.1 -lambda 1e6
export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}_T -J ${BASENAME}_T ./simulation_T -file ${FOLDER}/${BASENAME}_T -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 2.7 -ypos .9 -xpos ${P} -nu 0.000018382931330804 -sim falling -shape ellipse -semiAxisX 0.003125 -semiAxisY 0.025 -angle 0.2 -dlm ${D}
#export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}_T -J ${BASENAME}_T ./simulation_T -file ${FOLDER}/${BASENAME}_T -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 100. -rhoS 2.7 -ypos .9 -xpos ${P} -nu 0.0000065 -sim falling -shape ellipse -semiAxisX 0.0015625 -semiAxisY 0.0125 -angle 0.2 -dlm ${D}
#export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o ${BASENAME}_T_DLM10 -J ${BASENAME}_T_DLM10 ./simulation_T -file ${FOLDER}/${BASENAME}_T_DLM10 -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 2.7 -ypos .9 -xpos .15 -nu 0.000018382931330804 -sim falling -shape ellipse -semiAxisX 0.003125 -semiAxisY 0.025 -angle 0.2 -lambda 1e6 -dlm 10
# -tdump 0.1
#export OMP_NUM_THREADS=48;bsub -n 48 -W 40:00 -o ${BASENAME}${BPD}_F -J ${BASENAME}${BPD}_F ./simulation_F -file ${FOLDER}/${BASENAME}${BPD}_F -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.35 -ypos .9 -nu 0.0000190757887517147 -sim falling -shape ellipse -semiAxisX 0.003125 -semiAxisY 0.025 -angle 0.2 -tdump 0.1 -lambda 1e6
    cd -
done
done
done
done