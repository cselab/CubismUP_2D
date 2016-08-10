cd ../makefiles
make clean
make config=production poisson=split-fftw bc=box precision=double rk2=true particles=false dlm=true -j
cd ../launch

FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/GlowinskiLee_061215_Walls_box
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp Glowinski.sh ${FOLDER}

for D in 0.1
#1
#10
do
for C in 0.1
#0.01 0.05 0.1
do
for BPD in 16
#32
#4 8 16 32 64
do
	NAME=GlowinskiLee_Walls_bpdx${BPD}_CFL${C}_DLM${D}_box
	export OMP_NUM_THREADS=48
	bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file ${FOLDER}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*3] -radius 0.02 -tend 5. -rhoS 1.25 -ypos .67 -nu 0.00064 -sim falling -lambda ${D} -g 9.81
done
done
done