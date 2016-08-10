cd ../makefiles
make clean
make config=production poisson=split-fftw bc=box precision=double rk2=true particles=false dlm=true movingframe=true -j
cd ../launch

FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/GlowinskiLee_FixedLarge_151215_box
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp GlowinskiFixed.sh ${FOLDER}

for BPD in 8 16 32
#0.1 1
#10
do
for C in 0.05
#0.01 0.05 0.1
do
for D in 0.1
#4 8 16 32 64
do
	NAME=GlowinskiLee_bpdx${BPD}_CFL${C}_DLM${D}_box
	export OMP_NUM_THREADS=48
#bsub -n 48 -W 40:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file ${FOLDER}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*3] -radius 0.1 -tend 5. -rhoS 1.25 -ypos .2 -nu 0.00715541752799933 -sim falling -lambda ${D} -g 9.81
#bsub -n 48 -W 40:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file ${FOLDER}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*3] -radius 0.0208333333333333 -tend 5. -rhoS 1.25 -ypos .5 -nu 0.000680413817439772 -sim falling -lambda ${D} -g 9.8
#${FOLDER}/simulation -file ${FOLDER}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*3] -radius 0.02 -tend 5. -rhoS 1.25 -ypos .5 -nu 0.00064 -sim falling -lambda ${D} -g 9.81 -fdump 500
bsub -n 48 -W 40:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file ${FOLDER}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy ${BPD} -radius 0.0625 -tend 5. -rhoS 1.25 -ypos .5 -nu 0.00353553390593274 -sim falling -lambda ${D} -g 9.81
done
done
done