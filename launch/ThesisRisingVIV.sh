cd ../makefiles
make clean
make config=production poisson=split-fftw bc=box precision=double rk2=true particles=false dlm=true movingframe=true -j
cd ../launch

FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/RisingVIV_varNu_15022016
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp ThesisRisingVIV.sh ${FOLDER}

for BPD in 32 64
do
for C in 0.1
do
for D in 1
do
for M in 0.4 0.5 0.6
do
for N in 0.0001
do
NAME=RisingVIV_bpdx${BPD}_CFL${C}_M${M}_nu${N}_box
export OMP_NUM_THREADS=48
bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file ${FOLDER}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy ${BPD} -radius .025 -tend 50. -rhoS ${M} -ypos .8 -nu ${N} -sim falling -dlm ${D} -g 9.81
done
done
done
done
done